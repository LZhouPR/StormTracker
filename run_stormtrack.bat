@REM @echo off

call conda activate alex
cd /d E:/Simon/tracking/alexcrawford0927/final_version/program

rem Record the start time
set start_time=%time%

echo start time: %start_time%
rem Run the Python script
python C01_Reprojection.py
python C02_CycloneDetection.py
python C02_SystemDetection.py
python C03_Subset_Crossing_byGridMask_andLength.py
python C03_Subset_Crossing_byLatLon_andLength.py
python C04_CycloneStatSummary_AllStorms.py
python C05_Track_Aggregation_Append_events.py
python C05_Track_Aggregation_Append_intensity.py
python C05_Track_Aggregation_Append_rates.py 
python C05_Track_Aggregation_Append_trkden.py
python C06_ExportToCSV.py

echo "Completed all"

rem Record the end time
set end_time=%time%
echo end time: %end_time%

REM Calculate the total time elapsed
set /a hours=%end_time:~0,2%-%start_time:~0,2%
set /a minutes=%end_time:~3,2%-%start_time:~3,2%
set /a seconds=%end_time:~6,2%-%start_time:~6,2%

if %seconds% lss 0 (
    set /a minutes = %minutes% - 1
    set /a seconds = %seconds% + 60
)
if %minutes% lss 0 (
    set /a hours = %hours% - 1
    set /a minutes = %minutes% + 60
)
if %hours% lss 0 set /a hours = %hours% + 24

echo Total time elapsed: %hours% hours, %minutes% minutes, %seconds% seconds