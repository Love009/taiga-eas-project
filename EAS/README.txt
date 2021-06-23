0) "WriteCoordWithTimeErr.C"
    the function WriteCoordWithTimeErr() allows to
write txt-file of coordinates and sigmas for stations
(txt-file of coords is reqiered for this function, here it is "coordinatesGRANDE.txt").
This txt-file is necessary for the function Initialize() in the next C-files

1) "Recon.C"
    the function Recon(const char *FileName = FILE_TO_RECON) allows to 
reconstruct phi and theta angles of EASes and compare with the angles obtaned by HiSCORE. And draws (dt_exp - dt_theory) histograms,
where dt_theory is calculated using PhiRec and ThetaRec values.

2) "ToyMC.C"
    the function Generate() allows to
generate events and write them to the root-file "ToyMC.root"

3) "ReadGrande.C"
    the function ReadGrande(string fileToWrite) allows to
read raw data from Grande (address = DATA_DIRECTORY) and write to the root-file vector of events with >= N_ST_MIN stations per event
    the function main() allows to
read several directories and write to the root-file. One have to write dirictories adresses and the output file names.
! All the analysis results are saved in the folder "analysis_results"
    the function ReadBinary(string fileName, int stID, vector<Event> &evts) can be useful for debugging of one binary file.

4) "SelectGrHiscEvts.C"
    the function SelectGrHiscEvts() allows to
compare data from Grande (GRANDE_DATA_ADDRESS) and from HiSCORE (HiSCORE_DATA_ADDRESS) and write coinciding events to the root-file (OUTPUT_FILE_NAME)
! All analysis results are saved in the folder "analysis_results"

5) "ADCosc.C"
	the function DrawOsc_0_2() allows to
draw an oscillogram of one ADC-channel (0 or 1st) signal for the station using events in the data-file
similar is the function DrawOsc_1_3() for the 1st and 3d ADC-channels

6) "ADCpedestal.C"
	the function Pedestal(const char * FileName, int stID) allows to
drow distrebutions for the mean of 20 points ADC signal for 0 and 2 channels of the station, using events from the root-file
    the function Pedestals(const char * FileName = INPUT_FILE_NAME, const char * FileToWrite = OUTPUT_FILE_NAME) allows to
calculate the pedestal parameters for all the stations (0 and 2 channels) using the input file INPUT_FILE_NAME and write it to the txt-file OUTPUT_FILE_NAME

7) "DrawPedestal_vs_time.C"
    the function DrawPedestal_vs_time(int sepDrStNum = 1) allows to draw all the stations pedestal parameters dependence
    on one picture and for the station with number = sepDrStNum separately
! it uses all the files from the directory "analysis_results/pedestals" (pedestal261119.txt - an example of filename)
    
8) "ReadTreeOfEvts.h"
    the function ReadTreeOfEvts(const char *FileName, vector<Event> &evts) allows to
read a vector of Events from a root-file and put it in evts

9) "CenterMassMethod.C"
    the function FindEASCenter(const char *FileName = FILE_TO_READ) allows to
read the root-file and find the center of each EAS, compare the result with the coordinates obtaned by HiSCORE data analysis

10) "CheckTime.C"
    the function TimeInterval(const char * dir = DATA_DIRECTORY, , const char * hiscAdress = HiSCORE_DATA_ADDRESS) allows to
print to console HiSCORE and Grande operation time interval

11) "SignalSelect.C"
    the function CorrectedTree(const char * InitFileName = INITIALIZATION_FILE_NAME, const char * FileName = INPUT_FILE_NAME, const char * FileToWrite = OUTPUT_FILE_NAME, int threshold = 5, int lowLimit = -10) allows to
 write a tree from the file FileName with corrected time (using the start of ADC signal) to the file FileToWrite
*    The class Event based on the class Station is defined in the file "Tunka.h". Objects of the class Event are written in the all root-files in the trees.
For reading such type of the trees a dictionary is necessary. You can generate the dict using line 
//gInterpreter->GenerateDictionary("Event","Tunka.h")// in console. 

12) "FindTSystematic.C" 
    the function FindTSystematic(const char* fileName) alliwa to
find dt = (dt1,dt2, ..., dt19) systematic shifts for Grande stations, using data from the file fileName.

Also it is enough to include the file "ReadTreeOfEvts.h" in the script where the trees of events from the root-files are read.
** In the scrips the parameter N_STN_MIN presence - it is a minimal number of stations in event which is reqiered during the analysis.
Do not forget to change it in the head of a script depending on the aims.

*** PhiGen and ThetaGen data members (class Event, see "Tunka.h") are used either for MC simulated angles or for angles reconstructed from HiSCORE data.
______________________________________________________________________________________________________________
                              Stages of the analysis:
# Read data
ReadGrande.C
    main()
#Collect statistics for the trigger time corrections and throw away bad signals 
ADCpedestal.C
    Pedestals()
#Correct signals
SignalSelect.C
    CorrectedTree()
#Find maching for Grande and HiSCORE events
SelectGrHiscEvts.C
    SelectGrHiscEvts()
#Angle reconstruction
Recon.C
    Recon() (or one of main-s)
#EAS center reconstraction
CenterMassMethod.C
    FindEASCenter()
______________________________________________________________________________________________________________  
                                    EXAMPLE
                            for the 25st Dec 2019
                        in console using ROOT package

.L CheckTime.C
TimeInterval("/k2/DATA_GRANDE+ReX/2018-19/dec18/131218.rsg/", "/k1/prosin_group/HiSCORE/out19_0_500/out_1312.dat")

.L ReadGrande.C
ReadGrande("/k2/DATA_GRANDE+ReX/2018-19/nov18/010120.02.rsg/", "analysis_results/Gr010120.root");
.q
root -l

.L ADCpedestal.C
Pedestals("analysis_results/Gr270220.root", "analysis_results/pedestals/pedestal270220.txt")
.q
root -l

.L SignalSelect.C
CorrectedTree("analysis_results/pedestals/pedestal270220.txt", "analysis_results/Gr270220.root", "analysis_results/Gr270220Corrected.root" , 5, -10)
.q
root -l


.L SelectGrHiscEvts.C
SelectGrHiscEvts("analysis_results/Gr210220Corrected.root", "/k1/prosin_group/HiSCORE/out20_qfg_aug/out_0201.dat", "analysis_results/test.root")

.L Recon.C
Recon("analysis_results/Gr+Hisc301219.root")

//For several days it is possible to write:
//  Recon(FileNameDay1)
//  Recon(FileNameDay2)
//  Recon(FileNameDay3)
//The histograms will be joint for all the days. If you want to analise them separetely quit ROOT and run the script for each of them separately

.L CenterMassMethod.C
FindEASCenter("analysis_results/Gr+Hisc010120.root")
______________________________________________________________________________________________________________
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

root -l
ReadGrande("/k2/DATA_GRANDE+ReX/2019-20/dec19/210220.01.rsg/", "analysis_results/Gr291219.01.root");
ReadGrande("/k2/DATA_GRANDE+ReX/2019-20/dec19/291219.02.rsg/", "analysis_results/Gr291219.02.root");
ReadGrande("/k2/DATA_GRANDE+ReX/2019-20/dec19/291219.03.rsg/", "analysis_results/Gr291219.03.root");
ReadGrande("/k2/DATA_GRANDE+ReX/2019-20/dec19/291219.04.rsg/", "analysis_results/Gr291219.04.root");
JoinTrees("analysis_results/Gr291219.00.root", "analysis_results/Gr291219.01.root", "analysis_results/Gr291219.00-01.root");
JoinTrees("analysis_results/Gr291219.00-01.root", "analysis_results/Gr291219.02.root", "analysis_results/Gr291219.00-02.root");
JoinTrees("analysis_results/Gr291219.00-02.root", "analysis_results/Gr291219.03.root", "analysis_results/Gr291219.00-03.root");
JoinTrees("analysis_results/Gr291219.00-03.root", "analysis_results/Gr291219.04.root", "analysis_results/Gr291219.root");

JoinTrees("analysis_results/Gr221119.01.root", "analysis_results/Gr221119.02.root", "analysis_results/Gr221119.root");

.L ReadGrandeIrkutsk.C
ReadGrande("/k2/DATA_GRANDE+ReX/2019-20/feb20/170220.02.rsg/", "analysis_results/Gr170220.root");
ReadGrande("/k2/DATA_GRANDE+ReX/2019-20/feb20/270220.04.rsg/", "analysis_results/Gr270220.04.root");
JoinTrees("analysis_results/Gr270220.03.root", "analysis_results/Gr270220.04.root", "analysis_results/Gr270220.root");

.q
root -l

.L SelectGrHiscEvtsIrkutsk.C
SelectGrHiscEvts("analysis_results/Gr270220.root", "/k1/prosin_group/HiSCORE/out20_qfg_aug/out_2702.dat", "analysis_results/Gr+Hisc270220_newFolder.root")

.q
root -l
.L DrawSkyPoints.C
DrawSkyPoints("analysis_results/Gr+Hisc221119_newFolder.root", 0)
DrawSkyPoints("analysis_results/Gr+Hisc301219_newFolder.root", 0)
DrawSkyPoints("analysis_results/Gr+Hisc261119_newFolder.root", 0)
DrawSkyPoints("analysis_results/Gr+Hisc271119_newFolder.root", 0)
DrawSkyPoints("analysis_results/Gr+Hisc291119_newFolder.root", 0)
DrawSkyPoints("analysis_results/Gr+Hisc201219_newFolder.root", 0)
DrawSkyPoints("analysis_results/Gr+Hisc291219_newFolder.root", 0)
DrawSkyPoints("analysis_results/Gr+Hisc251219_newFolder.root", 0)
DrawSkyPoints("analysis_results/Gr+Hisc291219_newFolder.root", 0)
DrawSkyPoints("analysis_results/Gr+Hisc301219_newFolder.root", 0)
DrawSkyPoints("analysis_results/Gr+Hisc010120_newFolder.root")
______________________________________________________________________________________________________________
//========!!!!========
/k2/DATA_GRANDE+ReX/2019-20/oct19
/k1/prosin_group/HiSCORE/out_arxiv/out20_afg3_dec
//========!!!!========


.L ReadGrande.C
ReadGrande("/k2/DATA_T133+GRANDE/2018-19/jan19/060119/", "analysis_results/GrT133_060119.root")
.q
root -l

.L ADCpedestal.C
Pedestals("analysis_results/GrT133_060119.root", "analysis_results/pedestals/pedestal060119.txt")
.q
root -l

.L SignalSelect.C
CorrectedTree("analysis_results/pedestals/pedestal060119.txt", "analysis_results/GrT133_060119.root", "analysis_results/GrT133_060119Corrected.root" , 5, -10)
.q
root -l

.L Recon_old.C
Recon("analysis_results/GrT133_060119Corrected.root", 0)
.q
root -l

.L SelectGrHiscEvts.C
SelectGrHiscEvts("analysis_results/GrT133_060119.root", "analysis_results/out_0601.dat", "analysis_results/Gr+HiscT133_060119.root")

.q
root -l
______________________________________________________________________________________________________________
.L CheckTime.C
TimeInterval("/k2/DATA_GRANDE+ReX/2019-20/feb20/010220.rsg/", "/k1/prosin_group/HiSCORE/out20_qfg_aug/out_0102.dat")

TimeInterval("/k2/DATA_GRANDE+ReX/2019-20/feb20/030220.rsg/", "/k1/prosin_group/HiSCORE/out20_qfg_aug/out_0302.dat")

TimeInterval("/k2/DATA_GRANDE+ReX/2019-20/feb20/291219.01.rsg/", "/k1/prosin_group/HiSCORE/out20_qfg_aug/out_2102.dat")

TimeInterval("/k2/DATA_GRANDE+ReX/2019-20/feb20/291219.rsg/", "/k1/prosin_group/HiSCORE/out20_qfg_aug/out_2102.dat")
TimeInterval("/k2/DATA_GRANDE+ReX/2019-20/feb20/291219.01.rsg/", "/k1/prosin_group/HiSCORE/out20_qfg_aug/out_2102.dat")
TimeInterval("/k2/DATA_GRANDE+ReX/2019-20/feb20/291219.02.rsg/", "/k1/prosin_group/HiSCORE/out20_qfg_aug/out_2102.dat")

TimeInterval("/k2/DATA_GRANDE+ReX/2019-20/feb20/291219.rsg/", "/k1/prosin_group/HiSCORE/out20_qfg_aug/out_2102.dat")
TimeInterval("/k2/DATA_GRANDE+ReX/2019-20/feb20/291219.01.rsg/", "/k1/prosin_group/HiSCORE/out20_qfg_aug/out_2102.dat")

TimeInterval("/k2/DATA_GRANDE+ReX/2019-20/feb20/291219.rsg/", "/k1/prosin_group/HiSCORE/out20_qfg_aug/out_2102.dat")
TimeInterval("/k2/DATA_GRANDE+ReX/2019-20/feb20/291219.01.rsg/", "/k1/prosin_group/HiSCORE/out20_qfg_aug/out_2102.dat")
TimeInterval("/k2/DATA_GRANDE+ReX/2019-20/feb20/291219.02.rsg/", "/k1/prosin_group/HiSCORE/out20_qfg_aug/out_2102.dat")

TimeInterval("/k2/DATA_GRANDE+ReX/2019-20/feb20/291219.rsg/", "/k1/prosin_group/HiSCORE/out20_qfg_aug/out_2102.dat")
TimeInterval("/k2/DATA_GRANDE+ReX/2019-20/feb20/291219.01.rsg/", "/k1/prosin_group/HiSCORE/out20_qfg_aug/out_2102.dat")

TimeInterval("/k2/DATA_GRANDE+ReX/2019-20/feb20/210220.rsg/", "/k1/prosin_group/HiSCORE/out20_qfg_aug/out_2102.dat")
TimeInterval("/k2/DATA_GRANDE+ReX/2019-20/feb20/210220.01.rsg/", "/k1/prosin_group/HiSCORE/out20_qfg_aug/out_2102.dat")

TimeInterval("/k2/DATA_GRANDE+ReX/2019-20/feb20/250220.rsg/", "/k1/prosin_group/HiSCORE/out20_qfg_aug/out_2502.dat")
TimeInterval("/k2/DATA_GRANDE+ReX/2019-20/feb20/250220.01.rsg/", "/k1/prosin_group/HiSCORE/out20_qfg_aug/out_2502.dat")

TimeInterval("/k2/DATA_GRANDE+ReX/2019-20/feb20/270220.rsg/", "/k1/prosin_group/HiSCORE/out20_qfg_aug/out_2702.dat")
TimeInterval("/k2/DATA_GRANDE+ReX/2019-20/feb20/270220.01.rsg/", "/k1/prosin_group/HiSCORE/out20_qfg_aug/out_2702.dat")
TimeInterval("/k2/DATA_GRANDE+ReX/2019-20/feb20/270220.02.rsg/", "/k1/prosin_group/HiSCORE/out20_qfg_aug/out_2702.dat")
TimeInterval("/k2/DATA_GRANDE+ReX/2019-20/feb20/270220.03.rsg/", "/k1/prosin_group/HiSCORE/out20_qfg_aug/out_2702.dat")
TimeInterval("/k2/DATA_GRANDE+ReX/2019-20/feb20/270220.04.rsg/", "/k1/prosin_group/HiSCORE/out20_qfg_aug/out_2702.dat")
______________________________________________________________________________________________________________



