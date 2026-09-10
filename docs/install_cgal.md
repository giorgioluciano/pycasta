# Installazione di CGAL su Windows per pycasta

Questo documento descrive come installare CGAL e compilare il backend weighted Delaunay (`pocketlab`) per usare la triangolazione regular/weighted di CGAL invece del fallback SciPy.

## Prerequisiti

- **Python 3.11+** (es. Miniconda/conda)
- **Git**
- **Visual Studio Build Tools 2019** o VS 2019/2022 con supporto C++
- **CMake 3.18+**

## Passo 1: Clonare pocketlab

```bash
cd C:\github
git clone https://github.com/giorgioluciano/pocketlab.git
cd C:\github\pocketlab
```

## Passo 2: Installare CGAL e dipendenze con vcpkg

```bash
cd C:\github\pocketlab
git clone https://github.com/Microsoft/vcpkg.git
cd vcpkg
.\bootstrap-vcpkg.bat
.\vcpkg install cgal --triplet x64-windows
```

Questo installerà CGAL, Boost, GMP, MPFR e tutte le dipendenze in:
C:\github\pocketlab\vcpkg\installed\x64-windows

text

**Nota:** Questa operazione può richiedere 20-30 minuti e scaricare diversi GB.

## Passo 3: Compilare pocketlab

Torna in `C:\github\pocketlab` ed esegui:

```bash
cd C:\github\pocketlab
set CMAKE_TOOLCHAIN_FILE=C:\github\pocketlab\vcpkg\scripts\buildsystems\vcpkg.cmake
set VCPKG_TARGET_TRIPLET=x64-windows
pip install -e .
```

Se la build va a buon fine, vedrai:

```text
Successfully built pocketlab
Successfully installed pocketlab-0.1.0
```

## Passo 4: Copiare le DLL

Le DLL di runtime devono essere accessibili a Python. Copiale in `site-packages`:

```bash
copy C:\github\pocketlab\vcpkg\installed\x64-windows\bin\*.dll %CONDA_PREFIX%\Lib\site-packages\pocketlab\
```

Se usi un ambiente conda chiamato `pycasta`, il comando diventa:

```bash
copy C:\github\pocketlab\vcpkg\installed\x64-windows\bin\*.dll C:\Users\Giorgio\miniconda3\envs\pycasta\Lib\site-packages\pocketlab\
```

## Passo 5: Verificare

```bash
python -c "import pocketlab.core; print('OK:', pocketlab.__file__)"
```

Se stampa il percorso senza errori, l'installazione è riuscita.

## Passo 6: Usare in pycasta

Ora puoi eseguire script che usano `cgal_weighted_delaunay`:

```bash
cd C:\github\pycasta
python benchmark\compare_delaunay_backends.py
```

## Note importanti

1. **Non rimuovere la cartella `vcpkg`** dopo la build: le DLL servono a runtime.
2. **Se cambi ambiente Python**, devi ricopiare le DLL nel nuovo `site-packages`.
3. **Il backend weighted** è significativamente più lento di SciPy, ma produce una triangolazione topologicamente diversa che tiene conto dei raggi atomici.
4. **I power radii** non sono attualmente esposti dal wrapper Python; i raggi usati nei confronti sono circumsfere euclidee delle celle weighted.

## Troubleshooting

### ImportError: DLL load failed

Assicurati di aver copiato le DLL in `site-packages\pocketlab\`. Se il problema persiste, prova ad aggiungere temporaneamente al PATH:

```bash
set PATH=C:\github\pocketlab\vcpkg\installed\x64-windows\bin;%PATH%
python -c "import pocketlab.core; print('OK')"
```

### CMake non trova CGAL

Usa il toolchain file di vcpkg:

```bash
set CMAKE_TOOLCHAIN_FILE=C:\github\pocketlab\vcpkg\scripts\buildsystems\vcpkg.cmake
```

### Build fallisce con errori di compilazione

Verifica di avere Visual Studio Build Tools 2019 o successivo installato con il workload "Desktop development with C++".

