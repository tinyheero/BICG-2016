% Linux Software Installation Guide
% Andrew Roth

The following document gives a brief overview of the essential concepts required for installing software in Linux.

## Before We Begin
The following notes are about how to install software from the Linux command line. Most new versions of Linux include helpful graphical programs to install software. The steps below will typically only be necessary in the following cases.

1. You are working on a computer which does not have a GUI, for example a remote server.
2. You do not have administrative privileges on the computer, again most remote servers.
3. You are installing software which has not been nicely packaged for you Linux version. This last point applies to most Bioinformatics software.

## Variables
The first concept that needs to be covered is a system variable. When you are working in the shell it is possible to save certain commands and file paths in a shell variable.

The syntax for setting a variable is as follows.

```
export VAR=/some/path
```

It is very important that there be no space in between the = and the text. To see the value of the variable we use the `echo` command.

```
echo $VAR
```

> Note that whenever you reference the variable you need to put a `$` in front.

## System Variables

### `$PATH`
In order for Linux to run programs it needs to now which directories contain programs. This is specified by setting the `$PATH` variable. This variable stores a list of all directories which contain executable programs. The list is delimited by colons, `:`, and Linux will search the directories for programs in the order the directories are listed. 

For example setting

```
export PATH=/usr/bin:/usr/local/bin
```

would specify that only two directroies, `/usr/bin` and `/usr/local/bin` contain executable files. If you try to execute `myprogram` then Linux will first look in `/usr/bin` for the program, and then `/usr/local/bin`.

To see you $PATH variable execute the following.

```
echo $PATH
```

A useful way to add a directory to the current `$PATH` is as follows.

```
export PATH=/new/dir:$PATH
```

This will make it so that `/new/dir` is the first folder searched when you type a command. If the executable file can't be found in `/new/dir` the system will look in the other directories previously listed in `$PATH`.

You can also add new folders to the end of the `$PATH` so that they are the last folders searched.

```
export PATH=$PATH:/new/dir
```

### `$LD_LIBRARY_PATH`
In addition to knowing where the executable files are located, Linux also needs to know where to find the *library* files. Library files are program files, but they cannot be executed. They typically contain functions required by several programs, so they are placed in a single place and shared between programs.

The `$LD_LIBRARY_PATH` variable tells Linux where to find library files. It can be viewed and manipulated exactly like the `$PATH` variable.

### `$HOME`
The `$HOME` variable stores the path to your home directory. This can be useful when writing shell scripts.

### Setting System Variables
When you execute 

```
export PATH=/new/dir:$PATH
```

your `$PATH` variable will be updated. These updates will only last for the current session; once you log out your `$PATH` variable will be back to normal.

If you want to make your changes persist long term you can edit a special file to make Linux set you `$PATH` each time you login.

The two most common files are

- `~/.bashrc`
- `~/.bash_profile`

where `~` indicates your home folder. Which file you actually have and should edit depends on the version of Linux you use.

Adding the 

```
export PATH=/new/dir:$PATH
```
to these files will cause `/new/dir` to added to your `$PATH` whenever you log in.

## Installing Software

### Installation Using Make
Many programs for Linux are written in C or C++. The vast majority of these programs are installed in a similar way using a program called `make`.

The following steps should cover the installation of most programs which use `make`.

1. Download the file.
2. Decompress the file.
- If the file is `.tar.gz` use `tar -zxvf file_name.tar.gz`.
- If the file is `.tar.bz2` use `tar -jxvf file_name.tar.bz2`
3. Enter the newly created directory.
4. See if there are any files called `README` or `INSTALL` to read.
5. (not always) Execute `./configure`, this will check that the software can successfully be installed on the computer and configure any variables needed to compile the software.

> Note : To *compile* means to covert human readable source code, to machine executable code.

6. Execute `make`. This will turn the C/C++ code into an executable program.
7. (not always) Execute `make install`. This will copy any files to the system folders. Usually this means you can simply type `mypgrogram` to run the new software instead of giving the full path to where the program executable is located.

### Installing R Packages
Installing packages in `R` is relatively straightforward. The are two main mechanisms

1. Installation using CRAN and `install.packages`
2. Installation using Bionconductor and `bioclite`.

#### Installation Using CRAN
To install a package from CRAN, simply find the package name and execute `install.packages('package_name')` from `R`. Follow any prompts which come up.

#### Installation Using Bioconductor
To install packages from Bioconductor you simply execute something like the following from `R`.

```
source("http://bioconductor.org/biocLite.R")
biocLite("affy")
```

which would install the *affy* package. The first line loads the Bioconductor installation script, while the second line installs the package. Every Bioconductor package will include a snippet of code like this which you can copy and paste into your `R` environment.

### Installing Python Packages
Installing packages in Python is generally straightforward. The first step is usually to download the *source* file. One website which keeps a fairly extensive list of Python packages is <https://pypi.python.org/pypi>.

Once you have downloaded the package you will need to decompress the file and enter the newly created directory.

From this directory you can execute `python setup.py install` to install the software.

#### Additional Tools
One program which is quite useful is [pip](https://pypi.python.org/pypi/pip). If you install `pip` many Python packages can be installed simply using the following command. 	

```
pip install SomePackage
```

## Installing Software Without Administrative Privileges
Administrative privileges are essentially the ability to place files anywhere on the file system. This in turn means you can alter which programs are available to every user on the system.

Unless you are working on you personal computer, you are unlikely to to have administrative privileges on the computer. On servers, where many people use the same computer, it would be very dangerous for everyone to have administrative privileges. They could potentially delete other user's files, change the installed software without others knowing etc.. Usually a small set of trusted systems administration people will have administrative privileges (sometimes called *root* privileges for archaic Unix reasons).

Restricting admin rights to a small group of people is great from the perspective of maintaining a stable system. It can be a pain for users who want to install their own software though. Typically you need to get in touch with the systems admin and make a requests, which can be slow depending on your infrastructure.

Thankfully, Linux provides a powerful way to avoid this problem. The key idea is that users alter system variables such `$PATH` in their own accounts. This means they can have their own software stored in a place they can access, and it won't affect other users.

### Setting Up Personal Software
The general approach to installing software for your own use is to install the software somewhere you have permissions to write to. A common choice is to install software to `~/bin` and libraries to `~/lib` where `~` indicates your home directory. This approach is fine, but it will often lead to a lot of other folders such as `~/include` and `~/share` being created in your home folder.

I prefer to create a folder `~/install` and then place programs under `~/install/bin`. This means you will only have a single folder under your home directory. 

The approach I use is to create a folder `~/install/src` where I will download the *source* files, typically in `.tar.gz` or `.tar.bz2` format. For C/C++ programs the installation steps are very similar with one small modification. In step 5 instead of running

```
./configure
```

you should run

```
./configure --prefix=~/install
```

or sometimes the software dislikes using `~` so I specify the full path to my home folder

```
./configure --prefix=/home/andrew/install
```

where **andrew** is my user name.

The `--PREFIX=~/install` flag tells the configuration script to install all files for the software under `~/install`. So executable files will appear under `~/install/bin`, library files under `~/install/lib`, manual pages under `~/install/share` etc..

The rest of the procedure is the same as before.

To make the installed software available you need to change some system variables. The following lines can be place in your `.bashrc` or `.bash_profile` files under your home directory, and will make your personal software available when you login.

```
export INSTALL_DIR=$HOME/install

export PATH=$INSTALL_DIR/bin:$PATH

export LD_LIBRARY_PATH=$INSTALL_DIR/lib:$LD_LIBRARY_PATH.
```

This will define a variable `$INSTALL_DIR` which we use to set `$PATH` and `$LD_LIBRARY_PATH`.

#### Installing R Packages
Before you can install `R` packages, you need to install a custom version of `R`. `R` can be installed like any other C/C++ program as discussed above.

Once you have `R` installed you check which version the system in using by typing the following command.

```
which R
```

If everything is working correctly it would print something like

```
/home/andrew/install/bin/R
```

if it prints something like

```
/usr/local/bin/R
```

you are still using the system version, not your personal version.

Assuming you have installed a personal version of `R` you can then install packages just like before. The only difference is these packages will be installed for your personal version and not be visible to other users.

#### Installing Python Packages
To install Python packages you need to install a personal version of Python. Again this follows the same basic procedure for C/C++ programs outlined above.

Once Python is installed correctly, and you are using your personal version, then package installation is as described previously.

#### Considerations When Using Personal Software
Keeping your own version of software has many benefits. In particular you can avoid many problems if the system adminstrators update a package you rely on.

One problem which arises when using computer clusters is that your home directory on one computer, will typically not be visible on other computers in the cluster. As a result jobs run on the cluster will use the system version of the software instead of your personal versions.

The way to avoid this issue is to install your personal software somewhere all nodes in the cluster can see. There is usually a folder, let's call it `/cluster/data`, which is visible to all nodes. Such a folder has to exist otherwise nodes in the cluster can't share files. 

To install personal software there, I would create a folder `/cluster/data/andrew/install` and replace `~/install` with this path in the procedure outlined above.

#### Reproducible Research
> Warning: The following is a bit advanced and can be skipped.

When dealing with large projects which you plan to write papers about, it is essential you track the versions of software you use.

The simplest way to do this, is to write down the version number of every piece of software you use. This works well, but sometimes you might submit the paper and while it is out for review update some software. Now when you have to do new analyses to address reviewer comments you will have a different versions of the software installed. This can be annoying because it means redoing all the analyses again so everything has been analysed the same way and is comparable.

The procedure we outline above can easily be modified to keep *project* specific software. The basic idea is to create a project folder, say 

```
/cluster/data/andrew/projects/my_project
```

Under this folder you can create an `install` directory, say 

```
/cluster/data/andrew/projects/my_project/install
```

You can then alter you system variables to point to folders under this directory, for example

```
export $PATH=/cluster/data/andrew/projects/my_project/install/bin:$PATH
```

One issue is that you don't want to modify the `.bashrc` or `.bash_profile` files in your home directory to point to these paths, because you don't want to use the *project* software for day to day use. Maybe because the versions have fallen behind what you need for other projects. 

The way to side-step this is to create a file which sets the system paths, but is only called when you need to work on this project. 

For example I would create a file

```
/cluster/data/andrew/projects/my_project/install/bash_profile
```

which contains

```
export INSTALL_DIR=/cluster/data/andrew/projects/my_project/install

export PATH=$INSTALL_DIR/bin:$PATH

export LD_LIBRARY_PATH=$INSTALL_DIR/lib:$LD_LIBRARY_PATH
```

Now when I need to tell the system to use *project* software, I call the `source` command which reads the file.

```
source /cluster/data/andrew/projects/my_project/install/bash_profile
```

This will mean that I am using the *project* versions of any software I installed. When I log out though, I will be back to my *personal* versions.






