@echo off
for %%f in (unittests\src\*) do (
   cl.exe "%%f" /Iinclude /std:c++17 /EHsc /O2 /W4 && "%%~nf.exe" || exit /b
)
