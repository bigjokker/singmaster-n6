@echo on
REM i=7 Band II + Z-jump. No giant m.
REM Exact i=7 already did k<=200. This job: k=201..k_max except {K,K+1}.
REM N=F_16 F_17=1576239  K=F_14 F_17=602069  k_max=756136.
REM NOT next-prime from k through (k, N/2].
REM Pre-flight first; abort on failure.
REM Default 8 workers. Set FAMILY_WORKERS to change.
REM Checkpoint results\i7_sweep.jsonl; rerun this bat to resume.
REM json exists => refuse.
REM Do NOT start a second copy. Not i=8. Not E.

cd /d "%~dp0.."
if not exist results mkdir results

if exist results\i7_sweep.json (
  echo results\i7_sweep.json already exists. Not rerunning.
  goto :fail
)

echo.
echo === i=7 family sweep  Band II + Z-jump ===
echo     N=1576239  K=602069  k_max=756136
echo     Checkpoint: results\i7_sweep.jsonl
echo     Final:      results\i7_sweep.json
echo     Workers:    %FAMILY_WORKERS%  (default 8)
echo.

python scripts\family_sweep.py --i 7
if errorlevel 1 goto :fail

echo.
echo === i=7 sweep finished ===
pause
exit /b 0

:fail
echo i=7 sweep was not started, or JSON already exists, or pre-flight failed.
pause
exit /b 1
