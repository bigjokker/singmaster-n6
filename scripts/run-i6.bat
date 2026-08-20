@echo on
REM i=6 Band II + Z-jump. No giant m.
REM Exact i=6 already did k<=200. This job: k=201..k_max except {K,K+1}.
REM N=F_14 F_15=229970  K=F_12 F_15=87840  k_max=110318.
REM NOT next-prime from k through (k, N/2].
REM Pre-flight first; abort on failure.
REM Default 8 workers. Set FAMILY_WORKERS to change.
REM Checkpoint results\i6_sweep.jsonl; rerun this bat to resume.
REM json exists => refuse.
REM Do NOT start a second copy. Not i=7/8. Not E.

cd /d "%~dp0.."
if not exist results mkdir results

if exist results\i6_sweep.json (
  echo results\i6_sweep.json already exists. Not rerunning.
  goto :fail
)

echo.
echo === i=6 family sweep  Band II + Z-jump ===
echo     N=229970  K=87840  k_max=110318
echo     Checkpoint: results\i6_sweep.jsonl
echo     Final:      results\i6_sweep.json
echo     Workers:    %FAMILY_WORKERS%  (default 8)
echo.

python scripts\family_sweep.py --i 6
if errorlevel 1 goto :fail

echo.
echo === i=6 sweep finished ===
pause
exit /b 0

:fail
echo i=6 sweep was not started, or JSON already exists, or pre-flight failed.
pause
exit /b 1
