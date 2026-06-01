import re
import subprocess
from dataclasses import dataclass
from pathlib import Path

import h5py
import pytest
from astropy.io import fits


def fixtures_dir():
    return Path(__file__).resolve().parent / "fixtures"


@dataclass(frozen=True)
class RunCfg:
    exe: str
    opts: dict  # flag -> value (Path / str / None)
    template_file: Path

    @property
    def args(self):
        out = []
        for flag, val in self.opts.items():
            if val is None:
                continue
            out.extend([flag, str(val)])
        return out

    def missing_paths(self):
        missing = []
        for val in self.opts.values():
            if val is None:
                continue
            if not Path(val).exists():
                missing.append(str(val))
        if not Path(self.template_file).exists():
            missing.append(str(self.template_file))
        return missing

    @property
    def output_dir(self):
        return Path(self.opts.get("--output-dir"))


def fail_if_missing(paths, label="fixtures"):
    assert not paths, f"Missing required {label}: " + ", ".join(paths)


def skip_if_missing(paths, label="fixtures"):
    if paths:
        pytest.skip(f"Skipping because missing {label}: " + ", ".join(paths))


def run_cli(exe, args, cwd):
    cmd = [exe] + list(args)
    #print(cmd)
    return subprocess.run(cmd, capture_output=True, text=True, cwd=str(cwd))


def collect_schema(h5_path):
    """
    Template-driven HDF5 schema:
    path -> info dict (kind, dtype_kind, rank, shape, fields, field_specs)
    """
    schema = {}

    def dtype_kind(dt):
        if dt.names:
            return "compound"
        if dt.kind == "S":
            return "S"
        return dt.name

    with h5py.File(h5_path, "r") as f:

        def visitor(name, obj):
            path = name
            if isinstance(obj, h5py.Group):
                schema[path] = {"kind": "group"}
                return

            if isinstance(obj, h5py.Dataset):
                info = {
                    "kind": "dataset",
                    "dtype_kind": dtype_kind(obj.dtype),
                    "rank": len(obj.shape),
                    "shape": tuple(obj.shape),
                }

                if obj.dtype.names:
                    info["fields"] = list(obj.dtype.names)
                    field_specs = {}
                    for field in obj.dtype.names:
                        fdt = obj.dtype.fields[field][0]
                        base = fdt.base
                        fshape = tuple(fdt.shape)
                        if base.kind == "S":
                            bname = "S"
                        else:
                            bname = base.name
                        field_specs[field] = (bname, fshape)
                    info["field_specs"] = field_specs

                schema[path] = info

        f.visititems(visitor)

    return schema


def assert_schema_matches_template(produced_h5, template_h5):
    tmpl = collect_schema(template_h5)
    got = collect_schema(produced_h5)

    tmpl_paths = set(tmpl.keys())
    got_paths = set(got.keys())

    missing = sorted(tmpl_paths - got_paths)
    assert not missing, (
        f"Produced HDF5 missing paths from template (first 50): {missing[:50]}"
    )

    for path, tinfo in tmpl.items():
        if tinfo.get("kind") != "dataset":
            continue

        ginfo = got[path]
        assert ginfo.get("kind") == "dataset", (
            f"{path}: expected dataset, got {ginfo.get('kind')}"
        )

        assert tinfo["dtype_kind"] == ginfo["dtype_kind"], (
            f"{path}: dtype-kind mismatch (template {tinfo['dtype_kind']}, got {ginfo['dtype_kind']})"
        )

        assert tinfo["rank"] == ginfo["rank"], (
            f"{path}: rank mismatch (template {tinfo['rank']}, got {ginfo['rank']})"
        )

        assert tinfo.get("shape") == ginfo.get("shape"), (
            f"{path}: shape mismatch (template {tinfo.get('shape')}, got {ginfo.get('shape')})"
        )

        if tinfo["dtype_kind"] == "compound":
            assert tinfo.get("fields") == ginfo.get("fields"), (
                f"{path}: compound field list mismatch.\nTemplate: {tinfo.get('fields')}\nGot: {ginfo.get('fields')}"
            )
            assert tinfo.get("field_specs") == ginfo.get("field_specs"), (
                f"{path}: compound field dtype/shape mismatch.\n"
                f"Template: {tinfo.get('field_specs')}\nGot: {ginfo.get('field_specs')}"
            )


def parse_log(log_path, patterns):
    """
    patterns: dict name -> regex with one capture group (number)
    returns dict name -> int
    """
    text = log_path.read_text(errors="replace")
    assert text.strip(), "Log file is empty/whitespace"

    parsed = {}
    for key, pat in patterns.items():
        m = re.search(pat, text)
        assert m, f"Log missing required line for {key}: {pat}"
        g = m.group(1)
        parsed[key] = int(g) if g.isdigit() else 1
    return parsed


def collect_fits_schema(path):
    """
    Collect a lightweight "schema" of a FITS file:
    - HDU order + names
    - HDU type
    - cards count (header length)
    - for tables: nrows, ncols, total shape, and column (name, format, unit)
    - for images: data shape + dtype (if present)
    """
    path = Path(path)

    schema = {
        "path": str(path),
        "bytes": path.stat().st_size,
        "hdus": [],
    }

    with fits.open(path) as hdul:
        for idx, hdu in enumerate(hdul):
            h = {
                "index": idx,
                "name": hdu.name,
                "type": type(hdu).__name__,
                "header_cards": len(hdu.header),
            }

            # Tables
            if isinstance(hdu, (fits.BinTableHDU, fits.TableHDU)):
                nrows = 0 if hdu.data is None else len(hdu.data)
                ncols = 0 if hdu.columns is None else len(hdu.columns)

                h["table"] = {
                    "nrows": nrows,
                    "ncols": ncols,
                    "shape": (nrows, ncols),
                    "columns": [
                        {
                            "name": c.name,
                            "format": c.format,
                            "unit": c.unit,
                        }
                        for c in hdu.columns
                    ],
                }

            # Images (PRIMARY can be an image too)
            elif isinstance(hdu, (fits.ImageHDU, fits.PrimaryHDU)):
                if getattr(hdu, "data", None) is not None:
                    h["image"] = {
                        "shape": hdu.data.shape,
                        "dtype": str(hdu.data.dtype),
                    }
                else:
                    h["image"] = {
                        "shape": (),
                        "dtype": None,
                    }

            schema["hdus"].append(h)

    return schema


def assert_fits_schema_matches_template(produced_fits, template_fits):
    tmpl = collect_fits_schema(template_fits)
    got = collect_fits_schema(produced_fits)

    def _hdu_keys(schema):
        keys = []
        for k, v in schema.items():
            if k == "__file__":
                continue
            if isinstance(v, dict) and "kind" in v:
                keys.append(k)
        return set(keys)

    tmpl_keys = _hdu_keys(tmpl)
    got_keys = _hdu_keys(got)

    missing = sorted(tmpl_keys - got_keys)
    assert not missing, f"Produced FITS missing HDUs from template: {missing}"

    extra = sorted(got_keys - tmpl_keys)
    assert not extra, f"Produced FITS has extra HDUs vs template: {extra}"

    for key in sorted(tmpl_keys):
        tinfo = tmpl[key]
        ginfo = got[key]

        # HDU type must match
        assert tinfo["kind"] == ginfo["kind"], (
            f"{key}: HDU type mismatch (template {tinfo['kind']}, got {ginfo['kind']})"
        )

        # Data kind must match (table vs image)
        assert tinfo["data_kind"] == ginfo["data_kind"], (
            f"{key}: data_kind mismatch (template {tinfo['data_kind']}, got {ginfo['data_kind']})"
        )

        # Shapes must match
        assert tinfo.get("shape") == ginfo.get("shape"), (
            f"{key}: shape mismatch (template {tinfo.get('shape')}, got {ginfo.get('shape')})"
        )

        if tinfo["data_kind"] == "table":
            # Column count/rows already covered by shape, but keep messages clear
            assert tinfo["nrows"] == ginfo["nrows"], (
                f"{key}: nrows mismatch (template {tinfo['nrows']}, got {ginfo['nrows']})"
            )
            assert tinfo["ncols"] == ginfo["ncols"], (
                f"{key}: ncols mismatch (template {tinfo['ncols']}, got {ginfo['ncols']})"
            )

            # Columns: same names/order and same (format, unit)
            tcols = tinfo["columns"]
            gcols = ginfo["columns"]

            tnames = [c["name"] for c in tcols]
            gnames = [c["name"] for c in gcols]
            assert tnames == gnames, (
                f"{key}: column name/order mismatch.\nTemplate: {tnames}\nGot: {gnames}"
            )

            for tc, gc in zip(tcols, gcols, strict=True):
                cname = tc["name"]
                assert tc["format"] == gc["format"], (
                    f"{key}:{cname}: format mismatch (template {tc['format']}, got {gc['format']})"
                )
                assert tc["unit"] == gc["unit"], (
                    f"{key}:{cname}: unit mismatch (template {tc['unit']}, got {gc['unit']})"
                )
        else:
            # Image dtype comparison can be flaky for None/no-data; but safe for real images
            assert tinfo.get("dtype") == ginfo.get("dtype"), (
                f"{key}: dtype mismatch (template {tinfo.get('dtype')}, got {ginfo.get('dtype')})"
            )
