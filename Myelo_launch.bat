@echo off

REM cesta ke skriptum
set "WSL_DIST=Ubuntu-22.04"
set "WSL_SRC=/home/ciri/workflow/myelo/src"
set "WSL_MAIN=/home/ciri/work/myelo/src/project/xsvato01/archer_nf/main.nf"
set "WSL_CONFIG=/home/ciri/work/myelo/src/nextflow.config"
set "RUN_SH=RUN.sh"
set "NFX_BIN=/home/ciri/.local/bin/nextflow"

REM logy na disk D k runum
set "LOG_DIR_WIN=D:\ann\WORKFLOW\Myelo\runs"
set "LOG_DIR_WSL=/mnt/d/ann/WORKFLOW/Myelo/runs"

REM stahnout samplesheet
wsl.exe -d %WSL_DIST% bash -lc "cd %WSL_SRC% && wget -q -O samplesheet.csv 'https://docs.google.com/spreadsheets/d/1JieaHFcx7rXMdGSW5ex2kfrV4HdpVpofov6Kw2h5nto/export?format=csv' && cat samplesheet.csv"
echo ""
echo Samplesheet stazen

REM spustit workflow, log a reporty
echo Spousti se workflow
wsl.exe -d %WSL_DIST% bash -lc "set -euo pipefail; cd %WSL_SRC% && chmod +x %RUN_SH% && run_id=$(date +%%Y%%m%%d-%%H%%M%%S) && mkdir -p %LOG_DIR_WSL% && log=%LOG_DIR_WSL%/myelo_${run_id}.log && echo '>>> Log:' $log && ./%RUN_SH% -with-timeline %LOG_DIR_WSL%/timeline_\${run_id}.html -with-report %LOG_DIR_WSL%/report_\${run_id}.html -with-trace %LOG_DIR_WSL%/trace_\${run_id}.txt 2>&1 | tee \$log"
echo Hotovo dvacet!
echo Vysledky zde D:\ann\WORKFLOW\Myelo\runs
echo Cteni zde %LOG_DIR_WIN%\myelo_YYYYMMDD-HHMMSS.log
echo a zde %LOG_DIR_WIN%\report_*html, timeline_*.html, trace_*.txt
echo.
pause
