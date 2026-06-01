from sst1mpipe.io import get_pde_correction_factors

def test_load_pde_factors():

    pde_correction_factors = get_pde_correction_factors()
    keys = pde_correction_factors['mc_correction_for_PDE'].keys()
    assert 'tel_001' in keys
    assert 'tel_002' in keys
