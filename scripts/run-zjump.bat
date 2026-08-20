@echo on
REM Band I Z-jump remnant. NOT until-kill. NOT next-prime from k
REM through a Z-slab. Hang-guard p-k>20000 is obsolete.
REM Worklist: 89195 Stage-3 hang-guards + k=1000001..4126621.
REM Stragglers 4126622..4126646 already done. Family K, K+1 skipped.
REM Jump to first LIVE prime (NONE / PART-lower / p>N/2). Cap 12.
REM Same factorial kernel as Band II. numpy required.
REM Default 8 workers. Set ZJUMP_WORKERS to change.
REM Checkpoint jsonl after each chunk; rerun this bat to resume.
REM Hours, not overnight. Not Band II. Not E.
REM Do NOT start a second copy.
REM json exists => refuse.

cd /d "%~dp0.."
if not exist results mkdir results

if exist results\zjump.json (
  echo results\zjump.json already exists. Not rerunning.
  goto :fail
)

echo.
echo === Band I Z-jump  hang-guards + 1e6-K  cap 12 ===
echo     Spec:       docs\zjump-spec.md
echo     Checkpoint: results\zjump.jsonl
echo     Final:      results\zjump.json
echo     Workers:    %ZJUMP_WORKERS%  (default 8)
echo.

python scripts\zjump.py
if errorlevel 1 goto :fail

echo.
echo === Z-jump finished ===
pause
exit /b 0

:fail
echo Z-jump was not started, or JSON already exists, or pre-flight failed.
pause
exit /b 1
