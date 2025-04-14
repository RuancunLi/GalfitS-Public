# GalfitS

Galaxy imaging spectra fitting tool.

### Installation:

To install GalfitS, follow these steps:

1. Clone the repository:
   ```sh
   git clone https://github.com/RuancunLi/GalfitS-Public.git
   ```

2. Navigate to the cloned directory:
   ```sh
   cd GalfitS
   ```

3. Install the package and set up the alias and environment variable:
   ```sh
   python setup.py /path/to/your/data
   ```
   This command will install the package, add an alias (`galfits`), and set the `GS_DATA_PATH` environment variable in your shell configuration file (e.g., `.bashrc` or `.zshrc`). Restart your terminal or run `source ~/.bashrc` (or `~/.zshrc`) to apply changes.

4. Install the dependent packages:
   ```sh
   pip install -r requirements.txt
   ```

### Usage

You can run GalfitS with a configuration file like this :

```sh
galfits config.lyric --workplace /path/to/save/results
```

# Your Package Name

![Documentation](https://img.shields.io/badge/docs-online-brightgreen)
![GitHub Actions](https://github.com/ruancunli/GalfitS-Public/actions/workflows/docs.yml/badge.svg)

[📚 View Documentation](https://ruancunli.github.io/GalfitS-Public/)
