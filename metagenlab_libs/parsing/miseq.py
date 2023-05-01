import os
import re
import datetime
import pandas


def parse_sample_sheet(sheet_path):
    '''
    NextSeq 2023.01
    -----------------------------------------------------
    [Header]
    FileFormatVersion,2
    RunName,IMU-COV-082-WGS-059
    InstrumentPlatform,NextSeq1k2k
    IndexOrientation,Forward

    [Reads]
    Read1Cycles,151
    Read2Cycles,151
    Index1Cycles,8
    Index2Cycles,8

    [Sequencing_Settings]
    LibraryPrepKits,NexteraXT-CleanPlexFLEX

    [BCLConvert_Settings]
    SoftwareVersion,3.8.4
    AdapterRead1,CTGTCTCTTATACACATCT
    AdapterRead2,CTGTCTCTTATACACATCT
    FastqCompressionFormat,gzip
        
    
    MiSeq 2023.01
    -----------------------------------------------------
    [Header]
    IEMFileVersion,4
    Investigator Name,Sebastien Aeby
    Experiment Name,20200213_MTB
    Date,13.02.2020 OR Date,2021-09-02
    Workflow,GenerateFASTQ
    Application,FASTQ Only
    Assay,Nextera XT
    Description,
    Chemistry,Amplicon

    [Reads]
    151
    151

    [Settings]
    ReverseComplement,0
    Adapter,CTGTCTCTTATACACATCT

    [Data]
    Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
    193253,1052,,,N702,CGTACTAG,S504,AGAGTAGA,,
    191060,1053,,,N703,AGGCAGAA,S504,AGAGTAGA,,
    1801092178,1054,,,N705,GGACTCCT,S504,AGAGTAGA,,
    1907164214,1058,,,N710,CGAGGCTG,S504,AGAGTAGA,,
    20190718_positive_control,T+,,,N711,AAGAGGCA,S508,CTAAGCCT,,
    '''
    
    run2data = {}

    
    with open(sheet_path, 'r') as f:
        table_section = False
        reads_section = False
        run_date = None
        sample_table = []
        read_len = []
        for row in f:
            if row.startswith("Date"):
                run_date = row.rstrip().split(",")[1]
                run_date_split = re.split("\.|\/|-", run_date)
                if len(run_date_split[2]) == 4:
                    # 13.02.2020
                    run_date = f'{int(run_date_split[2])}-{int(run_date_split[1]):02d}-{int(run_date_split[0]):02d}'
                elif len(run_date_split[0]) == 4:
                    # 2020-02-13
                    run_date = f'{int(run_date_split[0])}-{int(run_date_split[1]):02d}-{int(run_date_split[2]):02d}'
            if row.startswith("Experiment Name") or row.startswith("RunName"):
                run_name = row.rstrip().split(",")[1]   
            if row.startswith("Library Prep Kit") or row.startswith("LibraryPrepKits") or row.startswith("Assay"):
                assay = row.rstrip().split(",")[1]
            if row.strip() == '[Reads]':
                reads_section = True
                continue
            if row.strip() in ['[Data]', '[Cloud_Data]']:
                print("TABLE!")
                table_section = True
                continue
            # if empty line or start of a new section
            if row.strip() == '' or row.strip()[0] == '[':
                table_section = False
                reads_section = False
            if table_section:
                sample_data = row.rstrip().split(",")
                sample_table.append(sample_data)
            if reads_section:
                if 'Read1Cycles' in row:
                    read_len.append(row.rstrip().split(",")[1])
                    
                elif 'Index1Cycles' in row:
                    continue
                else:
                    r_len = row.rstrip()
                    if r_len.isdigit():
                        read_len.append(r_len)
    if len(sample_data) > 0:
        df_samples = pandas.DataFrame(sample_table[1:], columns=sample_table[0])
    if not run_date:
        if 'ProjectName' in df_samples.columns:
            # IMU-COV-082-WGS-059_2022-12-16T12_18_14_e49d087
            run_date = re.match("[0-9]+-[0-9]+-[0-9]+",df_samples["ProjectName"][0].split(run_name)[1][1:]).group(0)
        else:
            run_date = None
    
    run2data["run_date"] = run_date
    run2data["run_name"] = run_name
    run2data["assay"] = assay
    if isinstance(read_len, list):
        read_len = read_len[0]
    run2data["NumCycles"] = read_len
    run2data["path"] = os.path.dirname(sheet_path)
    
    return run2data


def parse_runinfo_nextseq():
    '''
<?xml version="1.0"?>
<RunInfo Version="6">
	<Run Id="221216_VL00235_33_AACFTVGM5" Number="33">
		<Flowcell>AACFTVGM5</Flowcell>
		<Instrument>VL00235</Instrument>
		<Date>2022-12-16T12:49:26Z</Date>
		<Reads>
			<Read Number="1" NumCycles="151" IsIndexedRead="N" IsReverseComplement="N"/>
			<Read Number="2" NumCycles="8" IsIndexedRead="Y" IsReverseComplement="N"/>
			<Read Number="3" NumCycles="8" IsIndexedRead="Y" IsReverseComplement="Y"/>
			<Read Number="4" NumCycles="151" IsIndexedRead="N" IsReverseComplement="N"/>
		</Reads>
		<FlowcellLayout LaneCount="1" SurfaceCount="2" SwathCount="4" TileCount="4">
			<TileSet TileNamingConvention="FourDigit">
				<Tiles>
					<Tile>1_1101</Tile>
					<Tile>1_1102</Tile>
					<Tile>1_1103</Tile>
					<Tile>1_1104</Tile>
					<Tile>1_1201</Tile>
					<Tile>1_1202</Tile>
					<Tile>1_1203</Tile>
					<Tile>1_1204</Tile>
					<Tile>1_1301</Tile>
					<Tile>1_1302</Tile>
					<Tile>1_1303</Tile>
					<Tile>1_1304</Tile>
					<Tile>1_1401</Tile>
					<Tile>1_1402</Tile>
					<Tile>1_1403</Tile>
					<Tile>1_1404</Tile>
					<Tile>1_2101</Tile>
					<Tile>1_2102</Tile>
					<Tile>1_2103</Tile>
					<Tile>1_2104</Tile>
					<Tile>1_2201</Tile>
					<Tile>1_2202</Tile>
					<Tile>1_2203</Tile>
					<Tile>1_2204</Tile>
					<Tile>1_2301</Tile>
					<Tile>1_2302</Tile>
					<Tile>1_2303</Tile>
					<Tile>1_2304</Tile>
					<Tile>1_2401</Tile>
					<Tile>1_2402</Tile>
					<Tile>1_2403</Tile>
					<Tile>1_2404</Tile>
				</Tiles>
			</TileSet>
		</FlowcellLayout>
		<ImageDimensions Width="8208" Height="5541"/>
		<ImageChannels>
			<Name>green</Name>
			<Name>blue</Name>
		</ImageChannels>
	</Run>
</RunInfo>

    
    '''



def parse_runinfo(runinfo_path):
    '''
    MiSeq rininfo file
<?xml version="1.0"?>
<RunParameters xmlns:xsd="http://www.w3.org/2001/XMLSchema" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
  <EnableCloud>true</EnableCloud>
  <RunParametersVersion>MiSeq_1_1</RunParametersVersion>
  <CopyManifests>true</CopyManifests>
  <FlowcellRFIDTag>
    <SerialNumber>000000000-BY7WY</SerialNumber>
    <PartNumber>15028382</PartNumber>
    <ExpirationDate>2019-04-26T00:00:00</ExpirationDate>
  </FlowcellRFIDTag>
  <PR2BottleRFIDTag>
    <SerialNumber>MS6813532-00PR2</SerialNumber>
    <PartNumber>15041807</PartNumber>
    <ExpirationDate>2019-05-08T00:00:00</ExpirationDate>
  </PR2BottleRFIDTag>
  <ReagentKitRFIDTag>
    <SerialNumber>MS6925444-500V2</SerialNumber>
    <PartNumber>15033573</PartNumber>
    <ExpirationDate>2019-04-03T00:00:00</ExpirationDate>
  </ReagentKitRFIDTag>
  <Resumable>true</Resumable>
  <ManifestFiles />
  <AfterRunWashMethod>Post-Run Wash</AfterRunWashMethod>
  <Setup>
    <SupportMultipleSurfacesInUI>true</SupportMultipleSurfacesInUI>
    <ApplicationVersion>2.6.2.1</ApplicationVersion>
    <NumTilesPerSwath>14</NumTilesPerSwath>
    <NumSwaths>1</NumSwaths>
    <NumLanes>1</NumLanes>
    <ApplicationName>MiSeq Control Software</ApplicationName>
  </Setup>
  <RunID>180823_M03935_0099_000000000-BY7WY</RunID>
  <ScannerID>M03935</ScannerID>
  <RunNumber>99</RunNumber>
  <FPGAVersion>9.5.12</FPGAVersion>
  <MCSVersion>2.6.2.1</MCSVersion>
  <RTAVersion>1.18.54</RTAVersion>
  <Barcode>000000000-BY7WY</Barcode>
  <PR2BottleBarcode>MS6813532-00PR2</PR2BottleBarcode>
  <ReagentKitPartNumberEntered>15033573</ReagentKitPartNumberEntered>
  <ReagentKitVersion>Version2</ReagentKitVersion>
  <ReagentKitBarcode>MS6925444-500V2</ReagentKitBarcode>
  <PreviousPR2BottleBarcode />
  <PreviousReagentKitBarcode />
  <ExperimentName>20180823_divers</ExperimentName>
  <Chemistry>Amplicon</Chemistry>
  <Username>sbsuser</Username>
  <Workflow>
    <Analysis>GenerateFASTQ</Analysis>
  </Workflow>
  <EnableAnalysis>false</EnableAnalysis>
  <Reads>
    <RunInfoRead Number="1" NumCycles="251" IsIndexedRead="N" />
    <RunInfoRead Number="2" NumCycles="8" IsIndexedRead="Y" />
    <RunInfoRead Number="3" NumCycles="8" IsIndexedRead="Y" />
    <RunInfoRead Number="4" NumCycles="251" IsIndexedRead="N" />
  </Reads>
  <TempFolder>D:\Illumina\MiSeqTemp\180823_M03935_0099_000000000-BY7WY</TempFolder>
  <AnalysisFolder>D:\Illumina\MiSeqAnalysis\180823_M03935_0099_000000000-BY7WY</AnalysisFolder>
  <RunStartDate>180823</RunStartDate>
  <MostRecentWashType>PostRun</MostRecentWashType>
  <RecipeFolder>D:\Illumina\MiSeq Control Software\CustomRecipe</RecipeFolder>
  <ILMNOnlyRecipeFolder>C:\Illumina\MiSeq Control Software\Recipe</ILMNOnlyRecipeFolder>
  <SampleSheetName>MS6925444-500V2</SampleSheetName>
  <SampleSheetFolder>D:\Illumina\MiSeq Control Software\SampleSheets</SampleSheetFolder>
  <ManifestFolder>D:\Illumina\MiSeq Control Software\Manifests</ManifestFolder>
  <OutputFolder>D:\Illumina\MiSeqOutput\180823_M03935_0099_000000000-BY7WY</OutputFolder>
  <FocusMethod>AutoFocus</FocusMethod>
  <SurfaceToScan>Both</SurfaceToScan>
  <SaveFocusImages>false</SaveFocusImages>
  <SaveScanImages>true</SaveScanImages>
  <CloudUsername>sebastien.aeby@chuv.ch</CloudUsername>
  <RunManagementType>BaseSpaceCloud</RunManagementType>
  <CloudRunId>115694579</CloudRunId>
  <SendInstrumentHealthToILMN>true</SendInstrumentHealthToILMN>
    '''
    import xml.etree.ElementTree as ET
    root = ET.parse(runinfo_path).getroot()
    run_name = root.findall('ExperimentName')[0].text
    read_length = root.findall('Reads')[0].findall("RunInfoRead")[0].get("NumCycles")
    start_date = datetime.datetime.strptime( root.findall('RunStartDate')[0].text, '%y%m%d').strftime("%Y-%m-%d")

    run2data = {}
    run2data["samples"] = None
    run2data["run_date"] = start_date
    run2data["run_name"] = run_name
    run2data["assay"] = None
    run2data["reads"] = read_length
    
    return run2data

def format_date(datestr):
    from dateutil import parser as dateutilparser
    d = dateutilparser.parse(datestr)
    val = d.strftime("%Y-%m-%d")
    return val

def parse_runparemeters(runparameters_path):
    
    import xml.etree.ElementTree as ET
    
    root = ET.parse(runparameters_path).getroot()
    
    synonyms = {"FlowcellRFIDTag": "Flowcell", "Setup": ""}
    
    parameters = { 'NextSeq': ["InstrumentName",
                               "InstrumentType",
                               "InstrumentSerialNumber",
                               "FlowCellSerialNumber",
                               "FlowCellExpirationDate",
                               "FlowCellLotNumber", 
                               "FlowCellVersion",
                               "FlowCellMode", 
                               "RunElapsedTime", 
                               "ApplicationName"]
                  
                  , 'MiSeq': [("FlowcellRFIDTag", ["SerialNumber", "PartNumber", "ExpirationDate", "LotNumber"]), 
                              ("Setup", ["ApplicationName"]), 
                              "LibraryPrepKit",
                              "IndexKit"]
    }
    
    data = {}
    
    for machine in parameters:
        for parameter in parameters[machine]:
            if isinstance(parameter, str):
                if parameter in synonyms:
                    synonym = synonyms[parameter]
                else:
                    synonym = parameter
                try:
                    val = root.find(parameter).text
                    if "Date" in parameter:
                        val = format_date(val)
                    data[synonym] = val
                except AttributeError:
                    pass
            if isinstance(parameter, tuple):
                if parameter[0] in synonyms:
                    synonym = synonyms[parameter[0]]
                else:
                    synonym = parameter[0]  
                match = root.find(parameter[0])
                for subparam in parameter[1]:
                    print("subparam", subparam)
                    try:
                        val = match.find(subparam).text
                        if "Date" in subparam:
                            val = format_date(val)
                        data[f"{synonym}{subparam}"] = val
                    except:
                        pass

    # add MiSeq Manually
    if 'InstrumentType' not in data and len(data) > 0:
        data["InstrumentType"] = 'MiSeq'
    
    return data