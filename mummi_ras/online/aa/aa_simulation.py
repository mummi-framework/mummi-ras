from subprocess import PIPE, Popen
import os
from pprint import pprint
from time import sleep
import glob
import shutil

from mummi_ras.online.aa.aa_siminputs import AAinput
from mummi_core.utils.utilities import sys_call 

import MDAnalysis as mda

from logging import getLogger
LOGGER = getLogger(__name__)

class AAsimulation():

    def __init__(self, outpath, locpath):
        # Logistical variables
        self.popenenv = {
            "shell":                True,
            "universal_newlines":   True,
            "cwd":                  '.',
            #"env":                  os.environ,
        }

        self.is_local = True
        if outpath == locpath:
            self.is_local = False

        self._equi_process = None     # Initialize the process to None
        self._md_process = None  # Initialize the process to None
        self.trajname=""
        self.toponame="gro2amber.gro"
        self.trajtype="TRJ"
        self.outpath=outpath
        self.curIdx = "0"
        self.prev_filedate=0
        self.backup_count=0
        self.backup_rate=5
        # nstlim (number of step limit) of each simulation run
        self.nvt_nstlim=125000
        self.npt_nstlim=250000 #npt is separated into two runs
        # Set the max simtime of each sub-md sim very large - each sim will be stopped by the workflow
        self.md_nstlim=12500000  # 50 ns, 1500000 6 ns
        # max simulation time 12500000 * 4 fs = 50 ns
        self.max_simtime=12500000

    def setup(self):
        # TODO copy the required files to current working place
        LOGGER.info("Starting AMBER setup ....")
        #utils.sys_call("cp " + self.outpath + "/amber.* .")
        shutil.copyfile(self.outpath + "/amber.prmtop", "amber.prmtop")
        shutil.copyfile(self.outpath + "/amber.inpcrd", "amber.inpcrd")
        shutil.copyfile(self.outpath + "/gro2amber.gro", "gro2amber.gro")

        try:
            sys_call("cp " + self.outpath + "/md.*.rst .")
            sys_call("cp " + self.outpath + "/md.*.out .")
        except:
            LOGGER.info('No restart files are there. This is a new simulation')

    ## NVT
    def nvt(self):
        #pprint(os.environ)
        LOGGER.info("Starting NVT equilibrium ....")
        kwargs = self.popenenv

        amber_exe=shutil.which("pmemd.cuda")
        LOGGER.info("AMBER executable " + amber_exe )

        AAinput.nvt_input(self.nvt_nstlim)

        # pmemd.cuda -O -i nvt.in -p amber.prmtop -c amber.inpcrd -ref amber.inpcrd \
        # -r md.NVT.rst -x md.NVT.mdcrd -inf md.NVT.info -l md.NVT.log -o md.NVT.out
        cmd = [#"source /usr/gapps/mummi/lassen/amber18/amber.sh;",
               #"/usr/gapps/mummi/lassen/amber18/bin/pmemd.cuda -O",
               "pmemd.cuda -O", "-i nvt.in", "-p amber.prmtop", "-c amber.inpcrd",
                "-ref amber.inpcrd", "-r md.NVT.rst", "-x md.NVT.mdcrd", "-inf md.NVT.info",
                "-l md.NVT.log", "-o md.NVT.out"
               ]
        LOGGER.info("cmd = %s", " ".join(cmd))

        #self._equi_process = Popen(cmd, **kwargs)
        stdout = os.path.join(self.outpath, "nvt.stderr")
        stderr = os.path.join(self.outpath, "nvt.stdout")
        with open(stdout, "wb") as out, open(stderr, "wb") as err:
            self._equi_process = \
                Popen(' '.join(cmd), stdout=out, stderr=err, **kwargs)
            sleep(20)

            LOGGER.info("Process Running? %s", self._equi_process.poll())
            cuda_visible_devices = os.environ.get("CUDA_VISIBLE_DEVICES", None)
            LOGGER.info("CUDA_VISIBLE_DEVICES=%s", cuda_visible_devices)

            #wait until process complete
            self._equi_process.wait()

        if not os.path.isfile('md.NVT.rst'):
            raise Exception("md.NVT.rst is not there")

        # backup the calcualtion
        if self.is_local:
            try:
                sys_call("cp md.NVT.* aa_analysis_logger.log " + self.outpath)
            except:
                LOGGER.info('Backup files fail')

    ## NPT
    def npt(self, num):
        LOGGER.info("Starting NPT "+str(num)+" equilibrium ....")
        kwargs = self.popenenv

        amber_exe=shutil.which("pmemd.cuda")
        LOGGER.info("AMBER executable " + amber_exe )

        AAinput.npt_input(self.npt_nstlim)
        # pmemd.cuda -O -i npt.in -p amber.prmtop -c md.NVT.rst -ref md.NVT.rst \
        # -r md.NPT.rst -x md.NPT.mdcrd -inf md.NPT.info -l md.NPT.log -o md.NPT.out
        cmd =[]
        fname="md.NPT."+str(num)
        if num==1:
            cmd = [#"source /usr/gapps/mummi/lassen/amber18/amber.sh;",
                    "pmemd.cuda -O", "-i npt.in", "-p amber.prmtop", "-c md.NVT.rst",
                   "-ref md.NVT.rst", "-r "+fname+".rst", "-x "+fname+".mdcrd", "-inf "+fname+".info",
                   "-l "+fname+".log", "-o "+fname+".out"
                   ]
        elif num==2:
            cmd = [#"source /usr/gapps/mummi/lassen/amber18/amber.sh;",
                    "pmemd.cuda -O", "-i npt.in", "-p amber.prmtop", "-c md.NPT.1.rst",
                   "-ref md.NPT.1.rst", "-r md.0.rst", "-x "+fname+".mdcrd", "-inf "+fname+".info",
                   "-l "+fname+".log", "-o "+fname+".out"
                   ]

        LOGGER.info("cmd = %s", " ".join(cmd))
        #os.system(' '.join(cmd))

        stdout = os.path.join(self.outpath, "npt.stdout")
        stderr = os.path.join(self.outpath, "npt.stderr")
        with open(stdout, "wb") as out, open(stderr, "wb") as err:
            self._equi_process = \
                Popen(' '.join(cmd), stdout=out, stderr=err, **kwargs)

            sleep(5)
            LOGGER.info("Process Running? %s", self._equi_process.poll())
            cuda_visible_devices = os.environ.get("CUDA_VISIBLE_DEVICES", None)
            LOGGER.info("CUDA_VISIBLE_DEVICES=%s", cuda_visible_devices)

            #wait until process complete
            self._equi_process.wait()

        if num==1:
            if not os.path.isfile('md.NPT.1.rst'):
                raise Exception("md.NPT.1.rst is not there")
        if num==2:
            if not os.path.isfile('md.0.rst'):
                raise Exception("md.0.rst is not there")


        #utils.sys_call("cp md.NPT.rst md.0.rst")
        #shutil.copy("md.NPT.rst", "md.0.rst")
        # TODO backup
        if self.is_local:
            try:
                sys_call("cp md.NPT.* md.0.rst aa_analysis_logger.log " + self.outpath)
            except:
                LOGGER.info('Backup files fail')

    def _getRstFileIdx(self, rstfile):
        filenamelist=rstfile.split(".")
        if len(filenamelist) !=3 :
            raise Exception("Restart file "+rstfile+ " is not named correctly")
        return int(filenamelist[1])

    def _needEquilibrum(self, rstFile, outFile):
        if not os.path.isfile(rstFile):
            return True
        if not os.path.isfile(outFile):
            return True

        with open(outFile, 'r') as f:
            for line in f:
                if line[3:25]=='Final Performance Info':
                    return False

        return True

    def mdrun(self):
        # copy the require files to setup the calculations
        LOGGER.info("Start AMBER ....")
        if self.is_local:
            self.setup()

        # run the equilibrum - NVT, NPT1, and NPT2 if it is not run yet
        #if not os.path.isfile('md.NVT.rst'):
        if self._needEquilibrum('md.NVT.rst', 'md.NVT.out'):
            self.nvt()

        #if not os.path.isfile('md.NPT.1.rst'):
        if self._needEquilibrum('md.NPT.1.rst', 'md.NPT.1.out'):
            self.npt(1)

        #if not os.path.isfile('md.0.rst'):
        if self._needEquilibrum('md.0.rst', 'md.NPT.2.out'):
            self.npt(2)

        # restart file are named md.0.rst  md.1.rst  md.2.rst ...
        rstlist=glob.glob("md.[0-9]*.rst")
        LOGGER.info("MD restart file list {}".format(rstlist))
        idxList=[self._getRstFileIdx(item) for item in rstlist]
        idxList.sort()
        LOGGER.info("MD file index list {}".format(idxList))

        #find out the good rst file to restart MD
        while len(idxList)>0:
            rstfileidx = idxList.pop()
            rstfile="md."+str(rstfileidx)+".rst"
            if not os.path.isfile(rstfile):
                LOGGER.info("Restart file md."+str(rstfileidx)+".rst is not exsited")
            else:
                LOGGER.info("Use restart file md."+str(rstfileidx)+".rst for simulation")
                break
        else:
            raise Exception("idxList is empty after check rstart files")

        preIdx_int = rstfileidx
        preIdx=str(preIdx_int)
        curIdx=str(int(preIdx)+1)

        # these variables needed to be set for aa_analysis.
        # before the check of MD has reach the max time step in the previous run
        self.trajname="md."+curIdx+".mdcrd"
        self.trajnochamber="md."+curIdx+".nochamber.mdcrd"
        self.toponame="gro2amber.gro"
        self.trajtype="TRJ"
        self.curIdx=curIdx

        #check if MD has reach the max time step in the previous run
        if preIdx != "0":
            prv_simtime = self.get_simtime(preIdx)
            LOGGER.info("Prevous simulation time step - {}".format(prv_simtime))
            if prv_simtime > self.get_maxsimtime():
                LOGGER.info("Previous total time step - {} exceeds max value - {}, simulation STOP".format(prv_simtime,self.get_maxsimtime()))
                # if the current simulation is skip due to reaching max value then self.curIdx should be set to preDix
                self.curIdx = preIdx
                return

        # check the local trajectories
        if preIdx_int > 0:
            for i in range(preIdx_int):
                traj_local = "md."+str(i+1)+".mdcrd"
                traj_name=self.outpath+"/"+traj_local
                if os.path.isfile(traj_name):
                    if not os.path.isfile(traj_local): # skip if local already has traj
                        sys_call("ln -s "+traj_name+" .")
                else:
                    LOGGER.info("Previous trajectory file "+traj_name+" is not exsited")
                    raise Exception("Previous trajectory file "+traj_name+" is not exsited")
        #pprint(os.environ)
        LOGGER.info("mdrun Starting AMBER MD Production ....")

        ## MD production
        kwargs = self.popenenv
        AAinput.md_input(self.md_nstlim)
        amber_exe=shutil.which("pmemd.cuda")
        LOGGER.info("AMBER executable " + amber_exe )
        LOGGER.info("AMBER max simulation time steps " + str(self.get_maxsimtime()))
        # pmemd.cuda -O -i md.in -p amber.prmtop-c md.${PREV}.rst -r md.${THIS}.rst -x md.${THIS}.mdcrd -inf md.${THIS}.info -l md.${THIS}.log -o md.${THIS}.out
        cmd = [#"source /usr/gapps/mummi/lassen/amber18/amber.sh;",
                "pmemd.cuda -O",  "-i md.in", "-p amber.prmtop", "-c md."+preIdx+".rst",
                "-ref md."+preIdx+".rst", "-r md."+curIdx+".rst", "-x md."+curIdx+".mdcrd", "-inf md."+curIdx+".info",
                "-l md."+curIdx+".log", "-o md."+curIdx+".out"
               ]
        LOGGER.info("cmd = %s", " ".join(cmd))
        #os.system(' '.join(cmd))

        self._md_process = Popen(' '.join(cmd), **kwargs)

        sleep(5)
        LOGGER.info("Process Running? %s", self._md_process.poll())
        cuda_visible_devices = os.environ.get("CUDA_VISIBLE_DEVICES", None)
        LOGGER.info("CUDA_VISIBLE_DEVICES=%s", cuda_visible_devices)

        if self._md_process.poll():
            out, err = self._md_process.communicate()
            LOGGER.info("---------------- AMBER MD stdout --------------\n%s", out)
            LOGGER.error("---------------- AMBER MD stderr --------------\n%s", err)

    def backup(self):
        LOGGER.info("Backup the MD simulation data ....")
        '''
        #LOGGER.debug('checking for md_process: {}'.format(self._md_process.poll()))
        LOGGER.debug("AAsimulation.backup Process Running? %s", self._md_process.poll())
        while self._md_process.poll() is None:
            # sleep for 5 minutes
            sleep(300)
        '''

        try:
            sys_call("cp md."+self.curIdx+".* aa_analysis_logger.log " + self.outpath)
            #utils.sys_call("cp md."+self.curIdx+".??? aa_analysis_logger.log " + self.outpath)
            #utils.sys_call("rsync --append md." + self.curIdx + ".mdcrd " + self.outpath)
        except:
            LOGGER.info('Backup md.' + self.curIdx +'.* files fail')

    def mdcrd_checkstate(self):

        LOGGER.info("Check if " + self.trajname + " is ready for analysis")
        time_limit=0
        # check if  trajectory is written
        while not os.path.isfile(self.trajname):
            if time_limit>20:
                LOGGER.error("Trajectory {} is missing".format(self.trajname))
                raise Exception("Trajectory {} is missing".format(self.trajname))
            sleep(60)
            time_limit=time_limit+1

        if os.path.isfile(self.trajname):
            """
            # check if trajectory has frame by file size assume size of frame > 1 MB
            file_size = os.path.getsize(self.trajname)/1024.0/1024.0   # in MB
            time_limit=0
            while file_size < 1.0:
                if time_limit>20:
                    break
                sleep(60)
                file_size = os.path.getsize(self.trajname) / 1024.0 / 1024.0
                time_limit=time_limit+1
            """
            # check if trajectory is corrupted
            time_limit = 0
            while True:
                try:
                    aa = mda.Universe(self.toponame, self.trajname, format='NCDF', dt=100)
                    break
                except:
                    time_limit = time_limit + 1
                    if time_limit > 20:
                        LOGGER.error("Trajectory {} is corrupted".format(self.trajname))
                        raise Exception("Trajectory {} is corrupted".format(self.trajname))
                    sleep(60)

            # check if the trajectory has been updated
            curr_filedate = os.stat(self.trajname)[8]
            time_limit = 0
            while curr_filedate == self.prev_filedate:
                if time_limit > 20:
                    break
                sleep(60)
                time_limit = time_limit + 1
            self.prev_filedate = curr_filedate

            # backup file at backup_rate
            if self.is_local:
                if self.backup_count % self.backup_rate == 0:
                    self.backup()
                self.backup_count = self.backup_count + 1
            # end of backup

            # If reach the max sim time step, kill the simulation process
            cur_simtime=self.get_current_simtime()
            LOGGER.info("Current simulation time step - "+str(cur_simtime))
            if cur_simtime > self.get_maxsimtime():
                LOGGER.info("Current time step - "+str(cur_simtime)+" exceed max value - "+ str(self.get_maxsimtime()))
                self._md_process.kill()

            # Don't need nochamber
            '''
            inputs = """# cpptraj script to write as NoChamber
parm amber.prmtop
reference amber.inpcrd
trajin {}
trajout {}
            """.format(self.trajname, self.trajnochamber)

            with open('mdcrd_checkstate.in', 'w') as f:
                f.write(inputs)

            while time_limit < 20:
                try:
                    time_limit = time_limit + 1
                    utils.sys_call("cpptraj -p amber.prmtop -i mdcrd_checkstate.in > mdcrd_checkstate.log")

                    # backup file at backup_rate
                    if self.is_local:
                        if self.backup_count % self.backup_rate == 0:
                            self.backup()
                        self.backup_count = self.backup_count+1
                    # end of backup
                    break

                    # make sure
                    if self.get_simtime() > self.get_maxsimtime():
                         self._md_process.kill()

                except:
                    sleep(60)
            '''

        '''
        if not os.path.isfile(self.trajnochamber):
            raise Exception("AAsimulation.mdcrd_checkstate {} is not there".format(self.trajnochamber))
        '''

    def running(self):
        running = bool(
            self._md_process is not None
            and self._md_process.poll() is not None
            and self._md_process.poll() != 0
        )
        return running

    def stop(self):
        if self.running():
            self._md_process.kill()

    def get_trajname(self):
        return self.trajname

    def get_trajnochamber(self):
        return self.trajnochamber

    def get_toponame(self):
        return self.toponame

    def get_trajtype(self):
        return self.trajtype

    def get_process(self):
        return self._md_process

    def get_islocal(self):
        return self.is_local

    def get_current_simtime(self):
        cur_simtime=self.get_simtime(self.curIdx)
        return cur_simtime

    def get_simtime(self, idxStr):
        simtime = []
        with open("md."+idxStr+".out", 'r') as f:
            for line in f:
                if line[1:8] == "NSTEP =":
                    #print(line[30:44])
                    simtime.append(int(float(line[30:44]) * 250)-self.nvt_nstlim-self.npt_nstlim*2)
        if len(simtime) > 0:
            return simtime[-1]

        return 0

    def set_backuprate(self, rate):
        self.backup_rate=rate

    def set_maxsimtime(self, max_simtime):
        self.max_simtime=max_simtime

    def get_maxsimtime(self):
        return self.max_simtime
