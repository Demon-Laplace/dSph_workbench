# IniFile synchronization policy

Only models numbered `Fornax2047` and later are synchronized to the remote
repository. Models below `Fornax2047` must not be added in future updates.

The user-requested `IC_MODELS.csv` catalogue also includes parameter summaries
for older local models. This does not authorize uploading their INI or particle files.

## Model numbering

New models in the 50-prefixed sequence receive the next unused model number;
check the current repository before assigning it. The dual-beta candidate based
on the original Fornax5028 structure is **Fornax5029** (`IC_Fornax5029.ini`), with
inner beta=-0.2 and outer beta=+0.8. The next number is **Fornax5030**, followed by
Fornax5031 and onward. A new number identifies a parameter candidate and does
not imply that particles were generated or that a simulation was run.
