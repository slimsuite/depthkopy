#!/usr/bin/python

# See below for name and description
# Copyright (C) 2016 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
#  
# This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program; if not, write to 
# the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
#
# Author contact: <seqsuite@gmail.com> / School of Biotechnology and Biomolecular Sciences, UNSW, Sydney, Australia.
#
# To incorporate this module into your own programs, please see GNU Lesser General Public License disclaimer in rje.py

"""
Module:       DepthKopy
Description:  Single-copy read-depth and kmer based copy number analysis
Version:      0.2.1
Last Edit:    26/10/21
Copyright (C) 2021  Richard J. Edwards - See source code for GNU License Notice

Function:
    The function of this module will be added here.

    1. Check input files and programs.

    2. Check/generate BAM file. =bam

    2a. Fork out samtools depth/mpileup tmp files.
    2b. Compile and reshape temp files. =depfile

    3. Call DepthCopy.R to calculate CN and/or Regcheck. =regfile

    4. Check and report tables output.

Commandline:
    ### ~ Main DepthKopy run options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    seqin=FILE      : Input sequence assembly. (Must be *.fa or *.fasta for kmerself=T.) [None]
    basefile=FILE   : Root of output file names [diploidocus or $SEQIN basefile]
    scdepth=NUM     : Single copy ("diploid") read depth. If zero, will use SC BUSCO mode [0]
    bam=FILE        : BAM file of long reads mapped onto assembly [$BASEFILE.bam]
    reads=FILELIST  : List of fasta/fastq files containing reads. Wildcard allowed. Can be gzipped. []
    readtype=LIST   : List of ont/pb/hifi file types matching reads for minimap2 mapping [ont]
    dochtml=T/F     : Generate HTML DepthKopy documentation (*.docs.html) instead of main run [False]
    ### ~ Depth and Copy Number options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    busco=TSVFILE   : BUSCO full table [full_table_$BASEFILE.busco.tsv]
    quickdepth=T/F  : Whether to use samtools depth in place of mpileup (quicker but underestimates?) [False]
    depfile=FILE    : Precomputed depth file (*.fastdep or *.fastmp) to use [None]
    homfile=FILE    : Precomputed homology depth file (*.fasthom) to use [None]
    regfile=FILE    : File of SeqName, Start, End positions (or GFF) for read coverage checking [None]
    checkfields=LIST: Fields in checkpos file to give Locus, Start and End for checking [Hit,SbjStart,SbjEnd]
    gfftype=LIST    : Optional feature types to use if performing regcheck on GFF file (e.g. gene) ['gene']
    winsize=INT     : Generate additional window-based depth and CN profiles of INT bp (0 to switch off) [100000]
    winstep=NUM     : Generate window every NUM bp (or fraction of winsize=INT if <=1) [1]
    chromcheck=LIST : Output separate window analysis violin plots for listed sequences (or min size) + 'Other' []
    seqstats=T/F    : Whether to output CN and depth data for full sequences as well as BUSCO genes [True]
    ### ~ KAT kmer options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    kmerself=T/F        : Whether to perform additional assembly kmer analysis [True]
    kmerreads=FILELIST  : File of high quality reads for KAT kmer analysis []
    10xtrim=T/F         : Whether to trim 16bp 10x barcodes from Read 1 of Kmer Reads data for KAT analysis [False]
    ### ~ System options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    forks=X         : Number of parallel sequences to process at once [0]
    killforks=X     : Number of seconds of no activity before killing all remaining forks. [36000]
    forksleep=X     : Sleep time (seconds) between cycles of forking out more process [0]
    tmpdir=PATH     : Path for temporary output files during forking [./tmpdir/]
    ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import os, string, sys, time
slimsuitepath = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../')) + os.path.sep
sys.path.append(os.path.join(slimsuitepath,'libraries/'))
sys.path.append(os.path.join(slimsuitepath,'tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_db, rje_kat, rje_obj, rje_readcore, rje_rmd
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0.0 - Initial Compilation.
    # 0.1.0 - Added KAT kmers and self-homology. Renamed DepthKopy.
    # 0.1.1 - Fixed chromcheck setting bug.
    # 0.1.2 - Implemented kmerself=F toggle.
    # 0.2.0 - Added seqstats=T/F : Whether to output CN and depth data for full sequences as well as BUSCO genes [False]
    # 0.2.1 - Fixed occasional R BUSCO bug and renamed pngplots directory after basefile.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [ ] : Populate Module Docstring with basic info.
    # [Y] : Populate makeInfo() method with basic info.
    # [ ] : Add full description of program to module docstring.
    # [Y] : Create initial working version of program.
    # [N] : Add REST outputs to restSetup() and restOutputOrder()
    # [ ] : Add to SLiMSuite or SeqSuite.
    # [Y] : Add chromcheck=LIST : Output separate window analysis violin plots for listed sequences (or min size) + 'Other' []
    # [Y] : Add KAT read and self kmers and self-homology depth files that can be generated and added to depthcopy.R.
    # [ ] : Add check for KAT installation.
    # [ ] : Change pngplots/ to $BASEFILE.pngplots/
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copy_right) = ('DepthKopy', '0.2.1', 'October 2021', '2021')
    description = 'Single-copy read-depth based copy number analysis'
    author = 'Dr Richard J. Edwards.'
    comments = ['This program is still in development and has not been published.',rje_obj.zen()]
    return rje.Info(program,version,last_edit,description,author,time.time(),copy_right,comments)
#########################################################################################################################
def cmdHelp(info=None,out=None,cmd_list=[]):   ### Prints *.__doc__ and asks for more sys.argv commands
    '''Prints *.__doc__ and asks for more sys.argv commands.'''
    try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if not info: info = makeInfo()
        if not out: out = rje.Out()
        ### ~ [2] ~ Look for help commands and print options if found ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        cmd_help = cmd_list.count('help') + cmd_list.count('-help') + cmd_list.count('-h')
        if cmd_help > 0:
            rje.printf('\n\nHelp for {0} {1}: {2}\n'.format(info.program, info.version, time.asctime(time.localtime(info.start_time))))
            out.verbose(-1,4,text=__doc__)
            if rje.yesNo('Show general commandline options?',default='N'): out.verbose(-1,4,text=rje.__doc__)
            if rje.yesNo('Quit?'): sys.exit()           # Option to quit after help
            cmd_list += rje.inputCmds(out,cmd_list)     # Add extra commands interactively.
        elif out.stat['Interactive'] > 1: cmd_list += rje.inputCmds(out,cmd_list)    # Ask for more commands
        ### ~ [3] ~ Return commands ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        return cmd_list
    except SystemExit: sys.exit()
    except KeyboardInterrupt: sys.exit()
    except: rje.printf('Major Problem with cmdHelp()')
#########################################################################################################################
def setupProgram(): ### Basic Setup of Program when called from commandline.
    '''
    Basic Setup of Program when called from commandline:
    - Reads sys.argv and augments if appropriate
    - Makes Info, Out and Log objects
    - Returns [info,out,log,cmd_list]
    '''
    try:### ~ [1] ~ Initial Command Setup & Info ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        info = makeInfo()                                   # Sets up Info object with program details
        if len(sys.argv) == 2 and sys.argv[1] in ['version','-version','--version']: rje.printf(info.version); sys.exit(0)
        if len(sys.argv) == 2 and sys.argv[1] in ['details','-details','--details']: rje.printf('%s v%s' % (info.program,info.version)); sys.exit(0)
        if len(sys.argv) == 2 and sys.argv[1] in ['description','-description','--description']: rje.printf('%s: %s' % (info.program,info.description)); sys.exit(0)
        cmd_list = rje.getCmdList(sys.argv[1:],info=info)   # Reads arguments and load defaults from program.ini
        out = rje.Out(cmd_list=cmd_list)                    # Sets up Out object for controlling output to screen
        out.verbose(2,2,cmd_list,1)                         # Prints full commandlist if verbosity >= 2 
        out.printIntro(info)                                # Prints intro text using details from Info object
        cmd_list = cmdHelp(info,out,cmd_list)               # Shows commands (help) and/or adds commands from user
        log = rje.setLog(info,out,cmd_list)                 # Sets up Log object for controlling log file output
        return (info,out,log,cmd_list)                      # Returns objects for use in program
    except SystemExit: sys.exit()
    except KeyboardInterrupt: sys.exit()
    except: rje.printf('Problem during initial setup.'); raise
#########################################################################################################################
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: DepthKopy Class                                                                                               #
#########################################################################################################################
class DepthKopy(rje_readcore.ReadCore,rje_kat.KAT):
    '''
    DepthKopy Class. Author: Rich Edwards (2021).

    Str:str
    - BAM=FILE        : BAM file of reads mapped onto assembly [$BASEFILE.bam]
    - BUSCO=TSVFILE   : BUSCO full table [full_table_$BASEFILE.busco.tsv]
    - HomFile=FILE    : Precomputed homology depth file (*.fasthom) to use [None]
    - PAF=FILE        : PAF file of reads mapped onto assembly [$BASEFILE.paf]
    - RegFile=FILE    : File of SeqName, Start, End positions (or GFF) for read coverage checking [None]
    - SeqIn=FILE      : Input sequence assembly (sortnr/diphap modes) []
    - TmpDir=PATH     : Path for temporary output files during forking (not all modes) [./tmpdir/]

    Bool:boolean
    - KmerSelf=T/F    : Whether to perform additional assembly kmer analysis [True]
    - QuickDepth=T/F  : Whether to use samtools depth in place of mpileup (quicker but underestimates?) [False]

    Int:integer
    - DepAdjust=INT   : Advanced R density bandwidth adjustment parameter [8]
    - WinSize=INT     : Generate additional window-based depth and CN profiles of INT bp (0 to switch off) [100000]

    Num:float
    - WinStep=NUM     : Generate window every NUM bp (or fraction of winsize=INT if <=1) [1]

    File:file handles with matching str filenames

    List:list
    - CheckFields=LIST: Fields in checkpos file to give Locus, Start and End for checking [Hit,SbjStart,SbjEnd]
    - ChromCheck=LIST : Output separate window analysis violin plots for listed sequences (or min size) + 'Other' []
    - GFFType=LIST    : Optional feature types to use if performing regcheck on GFF file (e.g. gene) ['gene']
    - Reads=FILELIST  : List of fasta/fastq files containing reads. Wildcard allowed. Can be gzipped. []
    - ReadType=LIST   : List of ont/pb/hifi file types matching reads for minimap2 mapping [ont]

    Dict:dictionary    

    Obj:RJE_Objects
    - DB = rje_db.Database object
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.strlist = ['HomFile']
        self.boollist = ['KmerSelf','DocHTML']
        self.intlist = ['WinSize']
        self.numlist = ['WinStep']
        self.filelist = []
        self.listlist = ['ChromCheck']
        self.dictlist = []
        self.objlist = []
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True,setfile=True)
        self._setReadCoreAttributes()   # See rje_readcore
        self._setKatAttributes()        # See rje_kat
        self.setStr({})
        self.setBool({'KmerSelf':True,'SeqStats':True})
        self.setInt({'WinSize':100000})
        self.setNum({'WinStep':1})
        ### ~ Other Attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setForkAttributes()   # Delete if no forking
#########################################################################################################################
    def _cmdList(self):     ### Sets Attributes from commandline
        '''
        Sets attributes according to commandline parameters:
        - see .__doc__ or run with 'help' option
        '''
        for cmd in self.cmd_list:
            try:
                self._generalCmd(cmd)   ### General Options ### 
                self._forkCmd(cmd)  # Delete if no forking
                self._readCoreCmd(cmd)  # Will set all the core commands recognised.
                self._katCmd(cmd)       # Set kat commands recognised.
                ### Class Options (No need for arg if arg = att.lower()) ### 
                #self._cmdRead(cmd,type='str',att='Att',arg='Cmd')  # No need for arg if arg = att.lower()
                #self._cmdReadList(cmd,'str',['Att'])   # Normal strings
                #self._cmdReadList(cmd,'path',['Att'])  # String representing directory path 
                self._cmdReadList(cmd,'file',['HomFile'])  # String representing file path
                #self._cmdReadList(cmd,'date',['Att'])  # String representing date YYYY-MM-DD
                self._cmdReadList(cmd,'bool',['KmerSelf','DocHTML'])  # True/False Booleans
                self._cmdReadList(cmd,'int',['WinSize'])   # Integers
                self._cmdReadList(cmd,'float',['WinStep']) # Floats
                #self._cmdReadList(cmd,'min',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'max',['Att'])   # Integer value part of min,max command
                self._cmdReadList(cmd,'list',['ChromCheck'])  # List of strings (split on commas or file lines)
                #self._cmdReadList(cmd,'clist',['Att']) # Comma separated list as a *string* (self.str)
                #self._cmdReadList(cmd,'glist',['Att']) # List of files using wildcards and glob
                #self._cmdReadList(cmd,'cdict',['Att']) # Splits comma separated X:Y pairs into dictionary
                #self._cmdReadList(cmd,'cdictlist',['Att']) # As cdict but also enters keys into list
            except: self.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    ### <2> ### Main Class Backbone                                                                                     #
#########################################################################################################################
    def run(self):  ### Main run method
        '''
        # DepthKopy: Single-copy read-depth and kmer based copy number analysis

        DepthKopy is an updated version of the gregcnv methods of Diploidocus. Details to follow.
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('DocHTML'): return rje_rmd.docHTML(self)
            if not self.setup(): self.printLog('#ABORT','Problem during setup: aborted'); return False
            ### ~ [2] ~ Add main run code here ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.depthCopy(): return False
            return True
        except:
            self.errorLog(self.zen())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def setup(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.checkInput(busco=self.getNum('SCDepth')<=0): return False
            return True     # Setup successful
        except: self.errorLog('Problem during %s setup.' % self.prog()); return False  # Setup failed
#########################################################################################################################
    ### <3> ### Main Class Methods                                                                                      #
#########################################################################################################################
    def depthCopy(self):    ### Main Depth Copy method.
        '''
        Main Depth Copy method. This will generate the BAM if needed, generate a depth file, and then call the
        depthcopy.R script to perform the actual calculations.
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            basefile = self.baseFile()
            ## ~ [1a] Generate BAM and depth files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            depfile = self.getFastDep()         # This will generate the BAM file if needed
            if not rje.exists(depfile): raise IOError('Failed to create depth file')
            self.debug(depfile)
            ## ~ [1b] Generate homology BAM and depth files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            homfile = self.homFile()
            if homfile:
                self.printLog('#HOM', 'Self-minimap2 homology file: {0}'.format(homfile))
            ## ~ [1c] Generate KAT kmer files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.katKmers():
                self.setStr({'KatFile':'{0}.kat-counts.cvg'.format(basefile)})
            if self.getBool('KmerSelf') and self.selfKat():
                self.setStr({'KatSelf':'{0}.self.kat-counts.cvg'.format(basefile)})
            ### ~ [2] Perform Depth Copy analysis ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            depthcopy = self.callRscript()
            self.debug(depthcopy)

            #!# Add depthcopy.R vignettes of the different use cases:
            # - SCdepth
            # - RegCNV from GFF
            # - RegCNV from TSV
            return depthcopy
        except: self.errorLog('%s.depthCopy error' % self.prog())
#########################################################################################################################
    def homFile(self):  ### Checks and generates self-mapping homology file
        '''
        Checks and generates self-mapping homology file.
        << homfile if generated, else None
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getStrLC('HomFile') in ['f','false']:
                self.setStr({'HomFile':'None'})
                self.printLog('#HOM','No self-minimap2 homology analysis (homfile=F)')
                return None
            if not self.getStrLC('HomFile'):
                self.setStr({'HomFile':rje.baseFile(self.getStr('SeqIn'),strip_path=False)+'.fasthom'})
            fastdep = self.getStr('HomFile')
            ### ~ [2] Generate self BAM file if required ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            bamfile = rje.baseFile(self.getStr('SeqIn'),strip_path=False)+'.self.bam'
            makebam = self.force() or not rje.checkForFiles([fastdep],log=self.log)
            makebam = makebam and not rje.checkForFiles([bamfile],log=self.log)
            if makebam:
                imax = self.iMax()
                mapcmd = '{0} -aDP -t {1} -I {2}G -k19 -w19 -m200 {3} {3}'.format(self.minimap2(), self.threads(), imax, self.getStr('SeqIn'))
                samtools = self.samtools()
                mapcmd = mapcmd + ' | {0} view -b -h - | {0} sort - | tee {1} | {0} flagstat -'.format(samtools,bamfile)
                self.printLog('#SYS',mapcmd)
                for fline in os.popen(mapcmd).readlines():
                    if rje.chomp(fline) and fline[:1] != '0':
                        self.printLog('#SELF',rje.chomp(fline))
            if not self.checkBAMFile(bamfile, makeindex=True, bai=False, csi=False, needed=True): return None
            ### ~ [3] Generate self depth file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            quick = self.getBool('QuickDepth')
            self.setBool({'QuickDepth':True})
            fastdep = self.forkDepth(bamfile, fastdep, secondary=True, setstr=False)
            self.setBool({'QuickDepth': quick})
            return fastdep
        except: self.errorLog('%s.homFile error' % self.prog()); return None
#########################################################################################################################
    def makeOptionStr(self):  ### Returns the option string for the depthcopy.R Rscript.
        '''
        Returns the option string for the depthcopy.R Rscript.
        depfile=FILE [busco=FILE] [scdepth=INT] [regfile=FILE] [reghead=LIST] [gfftype=LIST] [winsize=INT] [winstep=NUM]
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            options = ['pngdir={0}.plots'.format(self.baseFile())]
            for cmd in ['depfile','busco','scdepth','regfile','gfftype','winsize','winstep','adjust','basefile',
                        'katfile', 'katself', 'homfile']:
                val =  self.getData(cmd)
                if val and val != 'None':
                    options.append('{0}={1}'.format(cmd,val))
            if self.list['ChromCheck']:
                options.append('chromcheck={0}'.format(','.join(self.list['ChromCheck'])))
            if self.debugging(): options.append('debug=TRUE')
            if self.getBool('SeqStats'): options.append('seqstats=TRUE')
            optionstr = ' '.join(options)
            return optionstr
        except: self.errorLog('%s.callRscript error' % self.prog())
#########################################################################################################################
    def callRscript(self,optionstr=''):  ### Calls the depthcopy.R Rscript and parses output.
        '''
        Calls the depthcopy.R Rscript and parses output.
        >> optionstr:str = String of options to give to depthcopy.R script.
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            rdir = '%slibraries/r/' % slimsuitepath
            if not optionstr: optionstr = self.makeOptionStr()
            ### ~ [2] Run ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            complete = False
            rcmd = 'Rscript {0}depthcopy.R {1}'.format(rdir, optionstr)
            self.printLog('#RCMD',rcmd)
            RCMD = os.popen(rcmd)
            rline = RCMD.readline()
            while rline:
                if '] #' in rline:
                    rline = ' '.join(rline.split('] ')[1:])
                    logstr = rje.chomp(rline).split()
                    self.printLog('{0}'.format(logstr[0].upper()), ' '.join(logstr[1:]))
                elif rline[:1] == '[':
                    self.verbose(v=0,text=rje.chomp(rline),newline=0)
                else:
                    self.verbose(v=1,text=rje.chomp(rline),newline=0)
                #!# Parse scdepth and other key points to printLog
                complete = complete or 'DepthCopy.R finished' in rline
                rline = RCMD.readline()
            RCMD.close()
            return complete
        except: self.errorLog('%s.callRscript error' % self.prog())
#########################################################################################################################
### End of SECTION II: DepthKopy Class                                                                                  #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION III: MODULE METHODS                                                                                         #
#########################################################################################################################

#########################################################################################################################
### END OF SECTION III                                                                                                  #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION IV: MAIN PROGRAM                                                                                            #
#########################################################################################################################
def runMain():
    ### ~ [1] ~ Basic Setup of Program  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    try: (info,out,mainlog,cmd_list) = setupProgram()
    except SystemExit: return  
    except: rje.printf('Unexpected error during program setup:', sys.exc_info()[0]); return
    
    ### ~ [2] ~ Rest of Functionality... ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    try: DepthKopy(mainlog,cmd_list).run()

    ### ~ [3] ~ End ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    except SystemExit: return  # Fork exit etc.
    except KeyboardInterrupt: mainlog.errorLog('User terminated.')
    except: mainlog.errorLog('Fatal error in main %s run.' % info.program)
    mainlog.endLog(info)
#########################################################################################################################
if __name__ == "__main__":      ### Call runMain 
    try: runMain()
    except: rje.printf('Cataclysmic run error: {0}'.format(sys.exc_info()[0]))
    sys.exit()
#########################################################################################################################
### END OF SECTION IV                                                                                                   #
#########################################################################################################################
