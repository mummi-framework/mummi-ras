# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
import os
import time
import glob
import MDAnalysis as mda
import shutil
import parmed as pmd
from argparse import ArgumentParser, RawTextHelpFormatter

import mummi_core
import mummi_core.utils.utilities as utils
from mummi_core.utils.naming import MuMMI_NamingUtils as Naming
from mummi_core.utils.utilities import add_timestamp
from mummi_core.utils.logger import init_logger

import mummi_ras
from mummi_ras.transformations.backmapping.grotop import GroTop
import mummi_ras.transformations.createsim.particlesim as particlesim

from logging import getLogger
LOGGER = getLogger(__name__)


class Backmapping (object):

    def __init__(self, args):
        LOGGER.info('Initialize ...')
        LOGGER.info('args: patch_id = {}, frame_id = {}'.format(args.patch_id, args.frame_id))

        self.patch_id = args.patch_id
        #self.aa2cg_patch_id = args.aa2cg_patch_id
        self.frame_id = Naming.cgframe(self.patch_id, args.frame_id)
        self.snapshot = Naming.snapshot(args.frame_id)
        self.fstype = args.fstype
        LOGGER.info('patch = {}, frame = {}, snapshot = {}'.format(self.patch_id, self.frame_id, self.snapshot))

        # read frame from
        LOGGER.info(Naming.dir_sim('cg', simname=self.patch_id))

        self.cg_dir=args.inpath
        if self.cg_dir == "":
            self.cg_dir = Naming.dir_sim('cg', simname=self.patch_id)

        self.res_dir = os.path.join(Naming.dir_res('backmapping'), "cg2aa_template_campaign3")

        self.outpath = args.outpath                    # created below
        if self.outpath == "":
            self.outpath = Naming.dir_sim('aa', simname=self.frame_id)

        # Also add full local patch if using local
        self.localpath = args.outlocal
        if self.localpath == "":
            # e.g.:         /tmp/var/backmapping
            # change to:    /tmp/var/backmapping/pfpatch_000000384620_f00000
            self.localpath = os.path.join(Naming.dir_local(), 'backmapping',
                                          add_timestamp(self.patch_id))
            #self.localpath = os.path.join(Naming.dir_local(), 'backmapping')
            #self.localpath = self.localpath + "/" + utils.add_timestamp(self.patch_id)

        LOGGER.info("final inpath   = {}".format(self.cg_dir))
        LOGGER.info("final outpath  = {}".format(self.outpath))
        LOGGER.info("final outlocal = {}".format(self.localpath))

        # check if gromacs is available
        gmx_exe = shutil.which("gmx")
        gmxd_exe = shutil.which("gmx_d")

        LOGGER.info("Single precision Gromacs executable: ({})".format(gmx_exe))
        LOGGER.info("Double precision Gromacs executable: ({})".format(gmxd_exe))
        if gmx_exe is None or gmxd_exe is None:
            raise OSError('Could not find gromacs executables')

        # create the outpath and local path at beginning
        os.makedirs(self.outpath, exist_ok=True)
        os.makedirs(self.localpath, exist_ok=True)

    def run(self):
        LOGGER.info('Run Backmapping ...')
        flag_success, flag_failure = Naming.status_flags('backmapping')

        gmx_exe=shutil.which("gmx")
        LOGGER.info("GROMACS single precision (Modified by C Neale) " + gmx_exe )

        gmxd_exe=shutil.which("gmx_d")
        LOGGER.info("GROMACS double precision (Modified by C Neale) " + gmxd_exe )

        try:
            iointerface = mummi_ras.get_io(self.fstype)

            begin = time.time()
            self.prepare()
            self.compute()
            backmap_end = time.time()
            LOGGER.info("Backmapping takes = {} seconds".format(backmap_end-begin))
            self.gromacs2amber()
            gro2amber_end=time.time()
            LOGGER.info("gromacs2amber takes = {} seconds".format(gro2amber_end - backmap_end))

            # don't need nochamber
            '''
            # get parmtop nochamber file for MDAnalysis
            self.ambernochamber()
            ambernochamber_end=time.time()
            LOGGER.info("ambernochamber takes = {} seconds".format(ambernochamber_end - gro2amber_end))
            '''

        except Exception as e:
            LOGGER.error(e)
            self.tarupall()
            fail_dir=utils.add_timestamp("failed")
            fail_dir_abs = "{}/{}".format(self.outpath, fail_dir)

            os.makedirs(fail_dir_abs, exist_ok=True)
            #utils.sys_call("mkdir -p {}/{}".format(self.outpath, fail_dir))
            shutil.copy("backmapping.tar.gz", fail_dir_abs)
            utils.sys_call("cp backmapping*.log {}/ || true".format(fail_dir_abs)) # save the python logs
            #utils.sys_call("cp backmapping.tar.gz {}/{}".format(self.outpath, fail_dir))
            # create an empty file in the outpath and call it “backmapping_failure”.
            #utils.sys_call("touch {}/backmapping_failure".format(self.outpath))

            iointerface.send_signal(self.outpath, flag_failure)
            raise e
        else:
            self.tarup()
            os.makedirs(self.outpath, exist_ok=True)
            #utils.sys_call("mkdir -p {}/".format(self.outpath))
            utils.sys_call("cp backmapping.tar.gz amber.prmtop amber.inpcrd gro2amber.gro backmapping*.log {} || true".format(self.outpath))
            # create an empty file in the outpath and call it “backmapping_success”.
            #utils.sys_call("touch {}/backmapping_success".format(self.outpath))
            iointerface.send_signal(self.outpath, flag_success)

        shutil.rmtree(self.localpath)

    def prepare(self):
        LOGGER.info('Backmapping preparation ...')
        #utils.sys_call("mkdir -p {}".format(self.localpath))
        os.makedirs(self.localpath, exist_ok=True)
        os.chdir(self.localpath)
        # copy files from template
        #print ('copy from ({})'.format(Naming.dir_res('backmapping')))
        utils.sys_call("cp -r {}/* .".format(self.res_dir))
        #print ('copy to ({})'.format(Naming.dir_local()))
        utils.sys_call("ln -s {}/topol.tpr my.input.tpr".format(self.cg_dir))
        #utils.sys_call("ln -s {}/system.top my.input.top".format(self.cg_dir))
        #utils.sys_call("sed 's/RAS_RAF/KRASCRAF/' {}/system.top > my.input.top".format(self.cg_dir))
        shutil.copyfile(self.cg_dir+"/system.top", "my.input.top")

        if not os.path.isfile('./my.input.tpr'):
            raise Exception("my.input.tpr is missing")
        if not os.path.isfile('./my.input.top'):
            raise Exception("my.input.top is missing")

        self.get_snapshot()

        # get reference pdbs from feedback-aa2cg
        # warning don't use full path, copy form current MUMMI_ROOT feedback-aa2cg dir)

        psim = particlesim.load_sim(Naming.dir_sim('cg', simname=self.patch_id))

        try:
            LOGGER.info('psim.prefRefPdbs list: {}'.format(psim.prefRefPdbs))
            feedback_aa_dir = Naming.dir_sim('feedback-aa')
            for file in psim.prefRefPdbs:
                pdb_file = os.path.join(feedback_aa_dir, os.path.split(file)[1])
                LOGGER.info('cp pdb feedback file: {}'.format(pdb_file))
                shutil.copy(pdb_file, ".")

                if "_CRD_" in pdb_file:
                    # check if HVR feedback and copy
                    #newPDBcopyList = glob.glob("{}/{}_HVR_*[0-9].pdb".format(feedback_aa_dir, pdb_file[0:15]))
                    newPDBcopyList = glob.glob("{}/{}*-HVR_*.pdb".format(feedback_aa_dir, pdb_file[0:15]))
                    LOGGER.info('Copy HVR feedback: {}'.format(newPDBcopyList))
                elif "_HVR_" in pdb_file:
                    # check if CRD feedback and copy
                    #newPDBcopyList = glob.glob("{}/{}_CRD_*[0-9].pdb".format(feedback_aa_dir, pdb_file[0:15]))
                    newPDBcopyList = glob.glob("{}/{}*-CRD_*.pdb".format(feedback_aa_dir, pdb_file[0:15]))
                    LOGGER.info('Copy CRD feedback: {}'.format(newPDBcopyList))
                if newPDBcopyList != []:
                    sortList = sorted(newPDBcopyList, key=lambda x: x[-19:-4])
                    newPDBcopy = sortList[-1][0:-4]  # note glob gives full file path
                    LOGGER.info('cp pdb feedback file for previus alt feedback: {}'.format(pdb_file))
                    LOGGER.info('Copy feedback List: {}'.format(newPDBcopy))
                    shutil.copy(newPDBcopy, ".")
        except AttributeError:
            LOGGER.info('psim.rasrafRefPdb list: {}'.format(psim.rasrafRefPdb))
            feedback_aa_dir = Naming.dir_sim('feedback-aa')
            for file in psim.rasrafRefPdb:
                pdb_file = os.path.join(feedback_aa_dir, os.path.split(file)[1])
                LOGGER.info('cp pdb feedback file: {}'.format(pdb_file))
                shutil.copy(pdb_file, ".")
        #cmd ="cp "+feedback_aa_dir +"/RAS_RAF_TERNARY_scfix-CYFpos-*.pdb ."
        #try:
        #    utils.sys_call(cmd)
        #except:
        #    LOGGER.info('No PDB reference files from feeback-aa are provided. ')

    def get_snapshot(self):
        LOGGER.info('get_snapshot from positions.tar ...')
        tarfile=self.cg_dir+"/positions.tar"
        if not os.path.isfile(tarfile):
            raise Exception("Backmapping.get_snapshot "+tarfile+" is missing")

        u=mda.Universe("my.input.tpr", tarfile, format='DDCMDTAR', snapshot=self.snapshot)

        for ts in u.trajectory:
            print(ts)
            u.atoms.write("my.input.pdb")

        if not os.path.isfile('./my.input.pdb'):
            raise Exception("Backmapping.get_snapshot my.input.pdb is missing")

    def compute(self):
        LOGGER.info('Compute backmapping ./dorun_campaign3_local.sh')
        utils.sys_call("./dorun_campaign3_local.sh")
        if not os.path.isfile('./final.top'):
            raise Exception("Backmapping.compute final.top is missing")
        if not os.path.isfile("./success.gro"):
            raise Exception("Backmapping.compute success.gro is missing")

    def gromacs2amber2fs(self):
        # for time step of 2 fs
        LOGGER.info('Backmapping.gromacs2amber2fs')
        gmx_top = pmd.load_file('gro2amber.top', xyz='gro2amber.gro')
        gmx_top.save('amber.prmtop', format='amber')
        gmx_top.save('amber.inpcrd', format='rst7')

        if not os.path.isfile('./amber.prmtop'):
            raise Exception("Backmapping.gromacs2amber4fs amber.prmtop is missing")
        if not os.path.isfile("./amber.inpcrd"):
            raise Exception("Backmapping.gromacs2amber4fs amber.inpcrd is missing")

    def gromacs2amber4fs(self):
        # for time step of 4 fs
        LOGGER.info('Backmapping.gromacs2amber4fs')
        # the converted amber inputs is amber.prmtop and amber.inpcrd
        utils.sys_call("./gro2amber.sh")

        if not os.path.isfile('./amber.prmtop'):
            raise Exception("Backmapping.gromacs2amber4fs amber.prmtop is missing")
        if not os.path.isfile("./amber.inpcrd"):
            raise Exception("Backmapping.gromacs2amber4fs amber.inpcrd is missing")

    def ambernochamber(self):
        input_str="""# cpptraj script to write as NoChamber
parm amber.prmtop
# parmstrip :SOL,NA,CL
parmwrite out amber.nochamber.prmtop nochamber
        """

        with open("parmtop_nochamber.in", 'w') as f:
            f.write(input_str)

        utils.sys_call("cpptraj -p amber.prmtop -i parmtop_nochamber.in > parmtop_nochamber.log")
        if not os.path.isfile('./amber.nochamber.prmtop'):
            raise Exception("Backmapping.ambernochamber amber.nochamber.prmtopis missing")

    def gromacs2amber(self):
        LOGGER.info('Backmapping.gromacs2amber')
        # fix GROMACS top and gro file before conversion
        grotop = GroTop("final.top")
        grotop.fixGroTop("gro2amber-int.top")
        # gro file fixing is in GroTop since it is related to the change of top file.
        grotop.fixGroFile("success.gro", "gro2amber-int.gro")
        # reorder zinc
        grotop.reorderZNtop("gro2amber-int.top", "gro2amber.top")
        grotop.reorderZNgro("gro2amber-int.gro", "gro2amber.gro")
        # end of fix GROMACS top and gro file

        if not os.path.isfile('./gro2amber.top'):
            raise Exception("Backmapping.gromacs2amber gro2amber.top is missing")
        if not os.path.isfile("./gro2amber.gro"):
            raise Exception("Backmapping.gromacs2amber gro2amber.gro is missing")

        # convert gromacs inputs to amber inputs using 4 fs timestep
        self.gromacs2amber4fs()


    def tarup(self):
        LOGGER.info('Backmapping.backup')
        files="my.input.tpr my.input.top my.input.pdb my.input.fix.pdb " \
              "topol.top backmapped.top 0-backward.gro " \
              "1-EMa.log 1-EMa.g96 output.mdrun.1-EMa " \
              "1-EMb.log 1-EMb.g96 output.mdrun.1-EMb " \
              "1-EMc.g96 1-EMc.log output.mdrun.1-EMc " \
              "2-EM.log 2-EM.gro output.mdrun.2-EM " \
              "3-MD-0.0002.gro 3-MD-0.0002.log output.mdrun.3-MD-0.0002 " \
              "3-MD-0.0005.gro 3-MD-0.0005.log output.mdrun.3-MD-0.0005 " \
              "3-MD-0.001.gro 3-MD-0.001.log output.mdrun.3-MD-0.001 " \
              "3-MD-0.002.gro 3-MD-0.002.log output.mdrun.3-MD-0.002  " \
              "aa.gro output.mdanalysis final.gro final.top " \
              "LANGEVIN.gro LANGEVIN.log output.mdrun.LANGEVIN LANGEVIN_pbcmol.gro " \
              "output.mdanalysis.3 success.gro backmapping*.log backmapping.out " \
              "output.restraint_satisfaction.3-MD-0.001 output.restraint_satisfaction.success " \
              "gro2amber.top gro2amber.gro empty.mdp full.top amber.prmtop amber.inpcrd tmp.tpr"

        if os.path.exists("backmapping.tar.gz"):
            os.remove("backmapping.tar.gz")
        #utils.sys_call("rm -f backmapping.tar.gz")

        try:
            utils.sys_call("tar -czf backmapping.tar.gz "+files+" >&/dev/null")
        except:
            LOGGER.info('tar -czf backmapping.tar.gz '+files)
            LOGGER.info('It is OK if the tar missing files')

    def tarupall(self):
        LOGGER.info('Backup all files for failed backmapping')

        if os.path.exists("backmapping.tar.gz"):
            os.remove("backmapping.tar.gz")

        try:
            utils.sys_call("tar -czf backmapping.tar.gz * >&/dev/null")
        except:
            LOGGER.info('tar -czf backmapping.tar.gz *')
            LOGGER.info('It is OK if the tar missing files')

        # modified from https://github.com/LLNL/maestrowf/blob/develop/maestrowf/maestro.py by Francesco Di Natale, dinatale3@llnl.gov.
def setup_argparser():
    """Set up the program's argument parser."""
    parser = ArgumentParser(prog="Backmapping", description="Backmapping: converts a Martini CG snapshot into an all-atom CHARMM \n"
                            "format and generates the necessary AMBER input files.",
                            formatter_class=RawTextHelpFormatter)

    parser.add_argument('--fstype', type=str, default="mummi", choices=["mummi", "simple","taridx"],
                        help='Chose how to write files/data (default: mummi)')

    parser.add_argument("-p", "--patch_id", type=str, required=True,
                        help="pfpatch: ID of patch selected by ML for backmapping")

    #parser.add_argument("-ap", "--aa2cg_patch_id", type=str, required=True,
    #                    help="pfpatch ID of patch form aa2cg feedback")

    parser.add_argument("-f", "--frame_id", type=int, required=True,
                    help="snapshot: ID of snapshot selected by ML for backmapping")

    parser.add_argument("-i", "--inpath", type=str, default="",
                        help="Input path, shoudl be a full path and if not set use env $MUMMI_ROOT/sims-cg/[sim name] (default in MuMMI)")

    parser.add_argument("-o", "--outpath", type=str, default="",
                        help="Output path, shoudl be a full path and if not set use env $MUMMI_ROOT/sims-aa/[sim name] (default in MuMMI)")

    parser.add_argument("-ol", "--outlocal", type=str, default="",
                        help="Local output path, shoudl be a full path and if set run everything in this path and only copy to --outpath in the end \n"
                        "WARNING, if --outpath not set (default in MuMMI) /sims-aa/[sim name]_[timestamp] will be added to this path.")

    parser.add_argument("-d", "--loglevel", type=int, default=2,
                        help="Level of logging messages to be output:\n"
                             "5 - Critical\n"
                             "4 - Error\n"
                             "3 - Warning\n"
                             "2 - Info (Default)\n"
                             "1 - Debug")

    parser.add_argument("-l", "--logpath", type=str, default="./",
                        help="Path for log information")

    parser.add_argument("-c", "--logstdout", action="store_true",
                        help="Log to stdout in addition to a file")

    parser.add_argument("--logfile", type=str, default="backmapping",
                        help="Name of backmapping logger")

    '''
    # comment for now
    parser.add_argument("-g", "--gromacs", type=str, default="/usr/local/gromacs/bin/",
                        help="GROMACS bin directory, used for e.g. grompp, mdrun minimization, make_ndx")

    parser.add_argument("-m", "--mpi", type=str, default="",
                        help="mpi run command used for all GROMACS MD mdrun's")

    parser.add_argument("-r", "--mdrunopt", type=str, default="",
                        help="Options to add to GROMACS mpi MD mdrun's")
    '''

    return parser


def main():

    # setup the parser
    # use those arguments to init backmapping
    # Set up the necessary base data structures to begin study set up.
    parser = setup_argparser()
    args = parser.parse_args()

    mummi_core.init()
    mummi_core.create_root()
    mummi_core.init_logger(argparser=args)

    b = Backmapping(args)
    b.run()

if __name__ == '__main__':
    main()
