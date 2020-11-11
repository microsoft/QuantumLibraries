PSI4_TEMPLATE = """
memory {N} GB

molecule {name} {
symmetry {symmetry}
{charge} {spin}
  {geometry}
}

set {
  basis {basis}
  scf_type {pk}
  reference {rhf}
  d_convergence {1e-8}
  e_convergence {1e-8}
}

If not asking broombridge:
{driver}('{method}')

If asking broombridge & driver=='energy'
e, wfn = {driver}('{method}', return_wfn=True)
fcidump(wfn, fname='fcidump', oe_ints=['EIGENVALUES'])
clean()
"""
