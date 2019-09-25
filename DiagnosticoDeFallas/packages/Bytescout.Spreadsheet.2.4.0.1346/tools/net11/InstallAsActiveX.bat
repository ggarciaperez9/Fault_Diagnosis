REM change current directory (required for Vista and higher)
@setlocal enableextensions
@cd /d "%~dp0"

REM coping Bytescout.Spreadsheet.dll into /System32/ as COM server libraries
copy Bytescout.Spreadsheet.dll %windir%\System32\Bytescout.Spreadsheet.dll
copy Bytescout.Spreadsheet.tlb %windir%\System32\Bytescout.Spreadsheet.tlb

REM register the dll as ActiveX library
%windir%\Microsoft.NET\Framework\v1.1.4322\regasm.exe %windir%\System32\Bytescout.Spreadsheet.dll /tlb:%windir%\System32\Bytescout.Spreadsheet.tlb /codebase
