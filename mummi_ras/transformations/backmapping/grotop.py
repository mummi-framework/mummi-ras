from collections import OrderedDict
import mummi_core.utils.utilities as utils

from logging import getLogger
LOGGER = getLogger(__name__)

class GroTop(object):

    def __init__(self, topfile):
        LOGGER.info('GroTop.init')
        self.numRmCL = 0
        self.contents=""
        with open(topfile, 'r') as f:
            for line in f:
                self.contents=self.contents+line
        # divide into sections
        sections = self.contents.split('[')
        if len(sections)<=1:
            raise Exception('Gromacs Top file {} does not have sections'.format(topfile))
        # Handle the header first
        self.sectionOD=OrderedDict()

        self.sectionOD['header']=sections[0]

        for secRaw in sections[1:]:
            sec = secRaw.split(']')
            if len(sec)!=2:
                raise Exception('Gromacs Top file section {} format is wrong'.format(sec))
            secName = sec[0].strip()
            self.sectionOD[secName]=sec[1]

    def getNumRmCL(self):
        return self.numRmCL

    def tofile(self, outfilename):
        LOGGER.info('GroTop.tofile')
        with open(outfilename, 'w') as f:
            for key, value in self.sectionOD.items():
                if key == 'header':
                    f.write(value)
                else:
                    f.write("[ "+key+" ]")
                    f.write(value)

    def molsection(self):
        LOGGER.info('GroTop.ismolsection')
        if 'molecules' not in self.sectionOD:
            raise Exception('Gromacs Top file molecules section is missing')
        sec = self.sectionOD['molecules']
        lines=sec.split('\n')
        return lines

    def rmzeromol(self):
        LOGGER.info('GroTop.rmzeromol')
        lines=self.molsection()
        newsec=""
        for line in lines:
            if not line.isspace():
                if len(line)>0 and line[0] != ';':
                    molvals=line.split()
                    if len(molvals)>1:
                        if int(molvals[1]) ==0:
                            continue
            newsec=newsec+line+"\n"

        self.sectionOD['molecules']=newsec

    def rmCLraf(self):
        LOGGER.info('GroTop.rmCLraf')
        rafList=['KRASCRAF', 'RAS_RAF' , 'RAS4A_RAF']
        # fine how many RASRAF or RAF molecules
        lines=self.molsection()
        numRaf=0
        for line in lines:
            if not line.isspace():
                if len(line) > 0 and line[0] != ';':
                    molvals=line.split()
                    if len(molvals)>1:
                        if molvals[0].strip() in rafList:
                            numRaf=numRaf+int(molvals[1])

        LOGGER.info('GroTop.rmCLraf - number of molecules contain RAF is {}'.format(numRaf))

        self.numRmCL=2*numRaf
        newsec=""
        for line in lines:
            if not line.isspace():
                if len(line) > 0 and line[0] != ';':
                    molvals=line.split()
                    if len(molvals)>1:
                        if molvals[0].strip() == "CL":
                            numCL=int(molvals[1])
                            numCL=numCL-self.numRmCL
                            if numCL < 0:
                                raise Exception('Number of CL is negative {}'.format(numCL))
                            newsec=newsec+"CL    {}\n".format(numCL)
                            continue
            newsec=newsec+line+"\n"

        self.sectionOD['molecules']=newsec

    def fixGroTop(self, outTopfile):
        LOGGER.info('GroTop.fixGroTop')
        self.rmzeromol()
        self.rmCLraf()
        self.tofile(outTopfile)

    # put the gro file fix here since it is related to the change of top file.
    def fixGroFile(self, inGrofile, outGrofile):
        LOGGER.info('GroTop.fixGroFile')
        numRmCL = self.numRmCL
        lines=[]
        with open(inGrofile, 'r') as f:
            for line in f:
                lines.append(line)

        if len(lines) <3:
            raise Exception('Gromacs gro file {} is wrong'.format(inGrofile))
        #fix the number of atom in Line[1] of gro file
        numAtom=int(lines[1].strip())-numRmCL
        lines[1]=str(numAtom)+"\n"

        #out put gro file
        with open(outGrofile, 'w') as of:
            of.write(lines[0])
            of.write(lines[1])
            for line in lines[2:-1]:
                if numRmCL > 0:
                    if len(line) < 15:
                        raise Exception('Gromacs gro file {} format is wrong: {}'.format(grofile, line))
                    resname = line[5:10].strip()
                    atmname = line[10:15].strip()
                    if resname == "CL" and atmname == "CL":
                        numRmCL = numRmCL -1
                        continue
                of.write(line)
            # write the last line
            of.write(lines[-1])

    def reorderZNtop(self, inTopfile, outTopfile):
        LOGGER.info('GroTop.reorderZNtop')
        # Use the reorder protein rascraf ITP file
        #utils.sys_call("sed -i -e 's/protein_rascraf.itp/neworder_protein_rascraf.itp/' {}".format(outTopfile))
        #utils.sys_call("sed -i -e 's/protein_rascraf.itp/protein_rascraf_reorder.itp/g;s/protein_ras4Acraf.itp/protein_ras4Acraf_reorder.itp/g' {}".format(outTopfile))
        substr1="protein_rascraf.itp"
        substr2="protein_ras4Acraf.itp"

        contents=[]
        with open(inTopfile, 'r') as f:
            for line in f:
                if substr1 in line:
                    s=line.replace(substr1, "protein_rascraf_reorder.itp")
                    contents.append(s)
                elif substr2 in line:
                    s=line.replace(substr2, "protein_ras4Acraf_reorder.itp")
                    contents.append(s)
                else:
                    contents.append(line)
        # write out use the same file name
        with open(outTopfile, 'w') as f:
            for line in contents:
                f.write(line)

    def reorderZNgro(self, inGrofile, outGrofile):
        LOGGER.info('GroTop.reorderZNgro')

        lines = []
        with open(inGrofile, 'r') as f:
            for line in f:
                lines.append(line)
        insertIdx=0
        count=1
        znList=[]
        for line in lines[2:-1]:
            count=count+1  # count start from 2
            if len(line) < 15:
                raise Exception('Gromacs gro file {} format is wrong: {}'.format(inGrofile, line))

            resname = line[5:10].strip()
            if resname == 'NMA':
                insertIdx = count + 1

            if resname == 'ZN2':
                delsize = len(znList)
                lines.pop(count-delsize)
                znList.append(line)

        for znVal in reversed(znList):
            lines.insert(insertIdx, znVal)

        with open(outGrofile, 'w') as of:
            for line in lines:
                of.write(line)
