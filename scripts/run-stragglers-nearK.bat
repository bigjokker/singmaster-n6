@echo on
REM 25 Band I stragglers k=4126622..4126646 (no prime in (k,K]).
REM Walk p > N/2 until each k is killed. One core. No giant m.
REM Expect a few minutes (first p=5401853, g~1.27e6, 25 columns).
REM Cap 40 primes so it cannot hang overnight.
REM Do not start a second copy. Do not start 1e6-K / Band II / E.

cd /d "%~dp0.."
if not exist results mkdir results

if exist results\stragglers_nearK.json (
  echo results\stragglers_nearK.json already exists. Not rerunning.
  goto :fail
)

echo.
echo === stragglers  i=8  k=4126622..4126646  p^>N/2 ===
echo.

python scripts\stragglers_nearK.py
if errorlevel 1 goto :fail

echo.
echo === stragglers finished ===
pause
exit /b 0

:fail
echo Stragglers was not started, or JSON already exists, or some k unkilled.
pause
exit /b 1
