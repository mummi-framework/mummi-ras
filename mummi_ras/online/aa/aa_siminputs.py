class AAinput:
    def __init__(self):
        pass

    @staticmethod
    def inputfile(name, irest, ntx, nstlim, temp0, tempi, ntb, ntp, pressure, restraint):
        # If ntwx or dt (Amber) are changed here, the dt (mda) = ntwx * dt (amber) in ps 
        # in the aa_analysis Universe load need to be changed.
        inputs="""A {simname} simulation
 &cntrl
    imin=0,          ! No minimization
    irest={simirest},         ! irest=1 for restart and irest=0 for new start
    ntx={simntx},           ! ntx=5 to use velocities from inpcrd and ntx=1 to not use them

    ! MD settings
    nstlim={simnstlim},
    dt=0.004,        ! time step (ps)

    ! Temperature control
    ntt=3,           ! Langevin dynamics
    temp0={simtemp0},     ! Target temperature
    tempi={simtempi},     ! Initial temperature
    gamma_ln=1.0,    ! Friction coefficient (ps^-1)
    ig=-1,

    ! Constant volume
    ntb={simntb},           ! 1=contant volume and 2=constant pressure
    ntp={simntp},           ! 1=isotropic, 2=anisotropic, 3=semi-isotropic w/ surften
    """.format(simname=name, simirest=irest, simntx=ntx, simnstlim=nstlim,
               simtemp0=temp0, simtempi=tempi, simntb=ntb, simntp=ntp)

        if pressure:
            inputs = inputs +"""barostat=2,      ! 1=Berendsen, 2=MC barostat
    pres0=1.01325,   ! Target external pressure, in bar
    taup=4,          ! Berendsen coupling constant (ps)
    comp=45,         ! compressibility

    ! Constant surface tension (needed for semi-isotropic scaling). Uncomment
    ! for this feature. csurften must be nonzero if ntp=3 above
    csurften=3,     ! Interfaces in 1=yz plane, 2=xz plane, 3=xy plane
    gamma_ten=0.0,  ! Surface tension (dyne/cm). 0 gives pure semi-iso scaling
    ninterface=2,   ! Number of interfaces (2 for bilayer)
            """

        inputs=inputs+"""
    ! SHAKE
    ntc=2,           ! Constrain bonds containing hydrogen
    ntf=2,           ! Do not calculate forces of bonds containing hydrogen

    ! Potential energy control
    cut=12.0,        ! nonbonded cutoff, in Angstroms
    fswitch=10.0,    ! for charmm.... note charmm-gui suggested cut=0.8 and no use of fswitch

    ! Control how often information is printed
    ntpr=25000,      ! Print energy frequency
    ntwx=25000,      ! Print coordinate frequency
    ntwr=25000,      ! Print restart file frequency
    ioutfm=1,        ! Write NetCDF format (always do this!)
    ntxo=2,          ! Write NetCDF format

    ! Wrap coordinates when printing them to the same unit cell
    iwrap=1,

    ! Set water atom/residue names for SETTLE recognition
    watnam='SOL',    ! Water residues are named TIP3
    owtnm='OW',      ! Water oxygens are named OH2
    hwtnm1='HW1',
    hwtnm2='HW2',
"""

        if restraint:
            inputs = inputs +"""
    ! position restraints
    ntr=1,
    restraint_wt=500.0,
    restraintmask='!(:SOL) & !(:NA) & !(:CL) & !(@H=)',
 """
        inputs = inputs +"""
 &end
 &ewald
    vdwmeth = 0,
 &end

"""
        return inputs

    @staticmethod
    def nvt_input(nstlim):
        #name, irest, ntx, nstlim, temp0, tempi, ntb, ntp, pressure, restraint
        inputs=AAinput.inputfile("NVT", 0, 1, nstlim, 310.0, 0.0, 1, 0, False, True)
        with open("nvt.in", 'w') as f:
            f.write(inputs)

    @staticmethod
    def npt_input(nstlim):
        inputs=AAinput.inputfile("NPT", 1, 5, nstlim, 310.0, 310.0, 2, 3, True, True)
        with open("npt.in", 'w') as f:
            f.write(inputs)

    @staticmethod
    def md_input(nstlim):
        inputs = AAinput.inputfile("MD", 1, 5, nstlim, 310.0, 310.0, 2, 3, True, False) # time step lim is changed to smaller one for debug
        with open("md.in", 'w') as f:
            f.write(inputs)

if __name__ == "__main__":
    AAinput.nvt_input()
    AAinput.npt_input()
    AAinput.md_input()
