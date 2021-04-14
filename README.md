Setup example:

`./ppj` is a project folder. The cpp code is stored in `./ppj/src/main.cpp` (you can change it's name)

Here is a minimal example

```bash
python3 -m venv dev         # create virtual env
source dev/bin/activate     # if you use bash alternatives
                            # use corresponding activate.* file if it exists
pip3 install ./ppj          # it now should install my_moudule 
python3 test.py             # should print ~5.9 and than fail because of type error
deactivate                  # to exit venv
```

venv (steps 1, 2, 5) are optional but highly recommended

