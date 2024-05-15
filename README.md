# TranStat
C code for implementing discrete-time chain-binomial transmission models
The configuration file "config.file" is where you specify all configuration parameters to run TranStat. 
For example datasets and configuration files, unzip the Examples.zip file. You will see four case studies.
Below are instructions on how to compile the C code

(1)	Mac/Linux

Assuming your mac OS or Linux has gcc facilities installed (majority should have), you can issue the following command:
gcc -O2 -lm -I /Users/yy70876/clib main.c -o main.exe
under the directory where the main source file “main.c” is located. This will produce an executable file, main.exe., which you will put (via file copy & paste) under the directory where your data and configuration files are.
Flag “-O2” is for optimization. You can change it to “-O1” or “-O3” for less or more optimization, respectively.  The flag “-I” tells the compiler about the “include” directory, i.e., where to find the needed head files (all the .h files). I put all the general-purpose head files under the directory “/Users/yy70876/clib”, so I need this flag. If you put the head files at the same directory as your “main.c”, then you don’t need this flag. Flag “-o” indicates the name of the output executable file. You can name it anything, not necessarily “main.exe”.
 
Tips: You don’t need to copy and paste “main.exe” each time. You can simply give the full path of “main.exe” under any folder, e.g., type “/Users/alex/transtat/main.exe” under folder “/Users/alex/data/”  

(2)	Windows

Go to https://www.bloodshed.net/ to download Dev-C++ 5. This is a great free C/C++ program editor and compiler. Installation is easy. If you are really not sure, check https://www.youtube.com/watch?v=NTkwZsUasXU. If you put head files in a different directory as I do, you can set the “Include” directory as follows:
1.	Go to the Tools menu.
2.	In the Tools menu, you should find and option called Compiler Options.
3.	On clicking Compiler Options, a dialog box should open with multiple tabs at the top which include Compiler (the default), settings, Directories and Programs. Click on Directories.
4.	Now you will find 4 new subtabs which include Binaries, Libraries, C includes, C++ includes. Click on “C includes”.
5.	Now you should see some buttons like Replace, Add etc. On top beside a text box you should see a button with the picture of a folder on it. Click on it.
6.	It should open another dialog box where you can choose the folder which contain the library files. After selecting the folder click Ok. Then click Add in the first dialog box. Then Ok again.

![image](https://github.com/yangyang-uf/TranStat/assets/25641021/5e9403eb-9da1-432a-8abd-68794f0c1664)

For more information, contact Prof. Yang Yang at yy70876@uga.edu
