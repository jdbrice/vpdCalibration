vpdCalibration
==============

VPD Calibration Project


##Workflow

1) Use the TOF Calibration nTuple Maker to produce nTuples in the TOF Calibration picoDst format. See TOF Calibration nTyple Maker for details.

2) Checkout this project into a clean working directory and build it:
```
$ cd some-working-directory
$ git clone <git-url>
$ cd vpdCalibration/bin
$ ./buildHere
```

3) Run the calibration with a valid configuration file 
```
$ ./vpd path/to/configuration/file.xml > path/to/log.txt &
```


## Configuration File

All configuration tags are case insensitive and must follow xml tag format
The default value is given for optional tags.

###Jobtype
* Default : calibration
  1. **calibrate**  
performs a geometry alignment job
  2. **paramReport**  
plots the parameter files given in the <paramInput> tag in the configuration file and compares them. Useful for comparing calibrations over time / different runs etc.
  3. **checkParams**
Readins in a parameter file then runs the calibration steps to produce qa plots.

###xVaraible
* Default : tof-tot
* The indep variable in the calibration procedure.
  1. **tof-tot** - Time of Flight's time-over-threshold ( 0 - 60 [ns] )
  2. **bbq-adc** - BBQ's adc value ( 0 - 4096 )
  3. **mxq-adc** - MXQ's adc value ( 0 - 4096 )

###yVaraible
* Default : tof-le
* The dep variable in the calibration procedure.
  1. **tof-le**  - Time of Flight's leading edge time [ns]
  2. **bbq-tdc** - BBQ's tdc value converted into ns 
  3. **mxq-tdc** - MXQ's tdc value converted to ns

###baseName
Default : ""
* The name to be prepended to all output files for easier record keeping. For instance give a basename of "run14AuAu14.6GeV" to make all root, PDF, etc appear as run14AuAu14.6GeV.{root, pdf, etc. } 

###rootOutput
* Default : "qa.root"
* The name specific to the root output file name. The full name will be baseName+rootOutput. The '.root' suffix will be added if needed.

###reportOutput
* Default : "qa.pdf"
* The name specific to the pdf report output file name. The full name will be baseName+reportOutput. The '.pdf' suffix should be specified.

###paramsOutput
* Default : "params.dat"
* The name specific to the data output file. The full name will be baseName+paramsOutput. The '.dat' suffix should be specified.

###paramsInput
* Default : {""}
* Vector of parameter file names to read in. Does not use the baseName prefix. If only one name is given and jobType=checkParams then the parameter file is applied and the calibration loop is started from there. If not or more than 1 filename is given, then the parameter files are plotted in overlay.

###paramsLegend
* Default : {""}
* Vector of parameter legend names to use when comparing parameter files. Does not use the baseName prefix.

###mapTriggerToTof
* Default : false 
* Determines wether or not the trigger (bbq and mxq ) channels will me mapped to tof channels. When producing final calibrations for the the trigger electronics, set this to false! 

###channelMap
* Default : ""
* The file containg the trigger to tof channel map. See the offline database for details on the format. If not given then <mapTriggerToTof> should be set to false.


###dataDir
* REQUIRED
* The full path to the directory containing the TOF calibration picoDsts

###maxFiles
* Default : 10000
* The maximum number of files to load from the <dataDir> directory for processing

###numIterations
* Default : 5
* The number of iterations to run the calibration procedure. Should be >= 4 for a good calibration. Usually use 8.

###minTOT
* Default : 10 [ns]
* The minimum tot value to consider in the slewing calibration.

###maxTOT
* Default : 40 [ns]
* The maximum tot value to consider in the slewing calibration.

###variableBinning
* Default : false 
* **True** - try to produce TOT bins that will provide equally distributed events in eahc TOT bin. Useful for lower statistics
* **False** - use fixed bin size detemined by the <numTOTBins> tag

###numTOTBins
* Default : 40
* provides the number of tot bins to use for variable or fixed binning. For fixed binning the bin size is simple the tot range / numTOTBins. Varialbe binning will produce numTOTBins but with varying sizes to accomidate the statistics in each tot region.

###removeOffset
* Default : true
* **True** - Calculates each detectors overall offset with respect to channel 1 on the west side, then removes it so all detector means are set to 0 relative to detector 1 on the west.
* **False** - Disallows the use of outlier rejection so <outlierRejection> should be false also.

###outlierRejection
* Default : true
* **True** - If the detector offsets are removed, the times can be used to calculate the VPD zVertex. This vertex is compared to the TPC zVertex and cuts are applied according the the <vzOutlierCut> tag.
* **False** - all times are used for calibration and no outlier rejection is performed.

###splineType
* Default : akima
* **akima** - use akima splines to fit the slewing curves and extract the corrections
* **cspline** - use a cubic slpine to fit the slewing curves and extract the corrections
* **linear** - use linear interpolation to fit the slewing curves and extract the corrections 
* **none** - use histogram bins to exctract the slewing corrections. Often causes discontinuities in the correction parameters.

###binMinPercent
* Default : 0.10 
* When using fixed binning, reject bins with too few events. threshold = (totalTotEvents/numTOTBins) * percent

###vzOutlierCut
* Default : { 40, 15, 8, 5} [cm]
* The vector of values corresponding to the vzCut for outlier rejection applied to the difference between the TPC zVertex and the calculated VPD zVertex. Each value corresponds to an iteration, first value is used on first iteration etc. The last value is used for all remaing iterations.

###avgNBackgroundCut
* Default : 10
* The threshold value applied to the avgN->FitSlicesY( function, firstBin, lastBin, avgNBackgroundCut )

###avgNTimingCut
* Default : { 2, 1, 0.6 } [ns]
* The timing cut applied when calculating the reduced average for the avgN detector resolution determination. Detector times that differe from the inclusive average by more than the cut will not be included in the reduced average.

###zeroCorrectionCut
* Default : 0.5 [ns]
* Used to remove wild fluctuations in the correction especially near the upper end of the tot range. After the <zeroStepN> iteration, if the differential correction between this step and the last is larger than <zeroCorrectionCut> then it is set to the value of teh previous bin with a "good" value.

###zeroStepN
* Default : 3
* The step at which the zeroCorrectionCut will begin being applied. Start counting at 0, should be atleast 2.

##Sample Configuration

An invalid configuration file will throw the error:
```
terminate called after throwing an instance of 'rapidxml::parse_error'
  what():  unexpected end of data
Abort
```
To fix, you should check for tags that aren't closed or other invalid xml.

```xml
<config>

	<!-- 	calibrate - for producing a new calibration 
			genReport - for reading in a calibration and producing qa -->
	<jobType>calibrate</jobType>

	<!-- output files -->
	<!-- The base name to use for all output of this job such as the root file, correction parameters etc.-->
	<baseName>a200</baseName>
	<!--	root, params, report filename appended to baseName -->
	<rootOutput>qa.root</rootOutput>
	<paramsOutput>params.dat</paramsOutput>
	<reportOutput>qa.pdf</reportOutput>

	<!-- Data input -->
	<!--
	<dataDir>/star/institutions/rice/jdb/run14/auau200/fastOffline/vpdCalibration/nTupler/dataDay75And80/output/</dataDir>-->
	<!--<dataDir>/star/institutions/rice/jdb/run13/pp510/tofCalibrationRun13/MuDstOutput/idealGeometry/output/</dataDir>-->
	<dataDir>/star/institutions/rice/jdb/run14/auau15/TofCalibration/geomAlign/t14Data0/output/</dataDir>
	<maxFiles>50</maxFiles>


	<!-- Calibration options -->
	<numIterations>4</numIterations>
	<variableBinning>false</variableBinning>
	<removeOffset>true</removeOffset>
	<outlierRejection>true</outlierRejection>
	<!-- 	linear: for linear interpolation
			cspline: for cubic spline interpolation
			akima:	for akima spline interpolation 
			anything else for bin based corrections -->
	<splineType>akima</splineType>
	<numTOTBins>40</numTOTBins>


	<vzOutlierCut>
		<v>40</v>
		<v>15</v>
		<v>8</v>
		<v>5</v> <!-- las value is for all remaining iterations -->
	</vzOutlierCut>

	<!-- Background cut applied to the avgN->FitSlicesY() when calculating detector resolution -->
	<avgNBackgroundCut>10</avgNBackgroundCut>
	<avgNTimingCut>
		<v>2</v>
		<v>1</v>
		<v>0.6</v> <!-- las value is for all remaining iterations -->
	</avgNTimingCut>

	<!-- After N steps if the difference in correction between this step and step N-1 is larger than zeroCorrectionCut, then the correction parameter is set to the value of the closest good value (by same criteria ) -->
	<zeroCorrectionCut>0.5</zeroCorrectionCut>
	<!-- The N above, the step at which the zeroCorrectionCut is applied. Should be greater than or equal to 1, start counting at 0, defaults to 2 -->
	<zeroStepN>3</zeroStepN>


</config>
```