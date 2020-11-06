The documentation is done in the thi - doc/ - folder. The documentation is written mainly in [mark down](https://www.markdownguide.org/) and uses [FORD](https://politicalphysicist.github.io/ford-fortran-documentation.html) on [GitHub](https://github.com/Fortran-FOSS-Programmers/ford) to combine MD and information created from the Fortran source files to create static HTML.

FORD can be installed via pip:

```
pip install ford
```

To generate the documentation simply do: 
```
ford getm_main.md
```

and then open your browser on ./html/index.html

