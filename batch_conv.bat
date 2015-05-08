ECHO

REM this creates a local variable for the convert function in IM
REM otherwise convert is confused with a window command
REM files in the directory are converted from pdf to tif

SETLOCAL EnableDelayedExpansion
SET IMCONV="C:\Program Files\ImageMagick-6.8.6-Q16\Convert"

FOR %%f IN (M:\docs\manuscripts\swmpr_manu\figure\*.pdf) do %IMCONV% -units PixelsPerInch -density 600 -resize 25%% "%%~f" "%%~dpnf.tif"
