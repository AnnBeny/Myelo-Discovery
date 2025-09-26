@echo off

set "WSL_DIST=Ubuntu-22.04"
set "WSL_SRC=/home/ciri/workflow/myelo/src"
set "RUN_SH=RUN.sh"

set "LOG_DIR_WIN=D:\ann\WORKFLOW\Myelo\runs"
set "LOG_DIR_WSL=/mnt/d/ann/WORKFLOW/Myelo/runs"

wsl.exe -d %WSL_DIST% bash -lc "set -e; cd %WSL_SRC% && chmod +x %RUN_SH% && run_id=$(date +%%Y%%m%%d-%%H%%M%%S) && log=%LOG_DIR_WSL%/myelo_${run_id}.log && echo '>>> Log:' $log && ./%RUN_SH% -with-timeline %LOG_DIR_WSL%/timeline_\${run_id}.html -with-report %LOG_DIR_WSL%/report_\${run_id}.html -with-trace %LOG_DIR_WSL%/trace_\${run_id}.txt 2>&1 | tee \$log"

pause
