@echo on
REM Band II p>N/2 image sweep. NOT a next-prime sweep from k.
REM Spec: docs\bandii-spec.md
REM Columns k=4126649..5182637 (1055989). Primes just above N/2.
REM Cap 14. Pre-flight §8 runs first and aborts on any failure.
REM Predictions: pass1 ~1.026e5, pass4 ~290, pass6 ~8.8, pass8 ~0.30.
REM numpy required. Default 8 workers. Set BANDII_WORKERS to change.
REM Checkpoint jsonl after each chunk; rerun this bat to resume.
REM Lunch break, not overnight. Not Band I. Not until-kill. Not E.
REM Do NOT start a second copy.
REM json exists => refuse.

cd /d "%~dp0.."
if not exist results mkdir results

if exist results\bandii_sweep.json (
  echo results\bandii_sweep.json already exists. Not rerunning.
  goto :fail
)

echo.
echo === Band II sweep  p^>N/2  1055989 columns  cap 14 ===
echo     Spec:       docs\bandii-spec.md
echo     Checkpoint: results\bandii_sweep.jsonl
echo     Final:      results\bandii_sweep.json
echo     Predict:    p1 ~1.026e5  p4 ~290  p6 ~8.8  p8 ~0.30
echo     Workers:    %BANDII_WORKERS%  (default 8)
echo.

python scripts\bandii_sweep.py
if errorlevel 1 goto :fail

echo.
echo === Band II sweep finished ===
pause
exit /b 0

:fail
echo Band II sweep was not started, or JSON already exists, or pre-flight failed.
pause
exit /b 1
