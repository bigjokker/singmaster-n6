@echo on
REM D only: Fibonacci i=8 extra-rep k<=300. One core, a few hours.
REM Start this after nearby_k333k-2M_de6.json looks clean.
REM See docs\campaign-log.txt (historical). Do not rerun unless you mean to.
REM
REM ALERT / NEW EXTRA REPRESENTATION = the result. Stop and paste it.
REM extras=[] / status=known_fibonacci is a certificate, not a failure.
REM Python still returns 0 on a discovery. JSON is written only when D ends.

cd /d "%~dp0"
if not exist results mkdir results

if not exist results\nearby_k333k-2M_de6.json (
  echo Missing results\nearby_k333k-2M_de6.json - nearby A is not done.
  goto :fail
)
if not exist results\fibonacci_i8.json (
  echo Missing results\fibonacci_i8.json - i=8 k^<=120 is not done.
  goto :fail
)
if exist results\fibonacci_i8_k300.json (
  echo results\fibonacci_i8_k300.json already exists. Not rerunning D.
  goto :fail
)

echo.
echo === D  Fibonacci i=8 extra-rep k^<=300 ===
echo     One core. Rebuilds the 3.1M-digit value, then inverts k=2..300.
echo     First it reprints C(10803704,4126647)=C(10803703,4126648).
echo     JSON appears only when this command ends. Do not start E.
echo.

python singmaster_intersect.py intersect --imin 8 --imax 8 --kextra 300 --json_out results\fibonacci_i8_k300.json
if errorlevel 1 goto :fail

echo.
echo === D finished ===
pause
exit /b 0

:fail
echo D was not started, or the search failed. Stopped here.
pause
exit /b 1
