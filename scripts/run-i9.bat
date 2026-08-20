@echo on
REM i=9 Band II + Z-jump. No giant m. Do NOT build m_9.
REM Modular already killed k=2..80. This job: k=81..k_max except {K,K+1}.
REM N=F_20 F_21=74049690  K=F_18 F_21=28284464  k_max=35522329.
REM Band II ~7.24e6 columns. Z-jump ~2.83e7 columns.
REM LARGER than i=8. Hours possible. F table ~300 MB per worker.
REM NOT next-prime from k through (k, N/2].
REM Pre-flight first; abort on failure.
REM Default 8 workers. Set FAMILY_WORKERS to change.
REM Checkpoint results\i9_sweep.jsonl; rerun this bat to resume.
REM json exists => refuse.
REM Do NOT start a second copy. Not i=6/7/8. Not E.

cd /d "%~dp0.."
if not exist results mkdir results

if exist results\i9_sweep.json (
  echo results\i9_sweep.json already exists. Not rerunning.
  goto :fail
)

echo.
echo === i=9 family sweep  Band II + Z-jump ===
echo     N=74049690  K=28284464  k_max=35522329
echo     This is LARGER than i=8. Hours possible.
echo     Checkpoint: results\i9_sweep.jsonl
echo     Final:      results\i9_sweep.json
echo     Workers:    %FAMILY_WORKERS%  (default 8)
echo.

python scripts\family_sweep.py --i 9
if errorlevel 1 goto :fail

echo.
echo === i=9 sweep finished ===
pause
exit /b 0

:fail
echo i=9 sweep was not started, or JSON already exists, or pre-flight failed.
pause
exit /b 1
