
External program usage
----------------------


Requirements for external program

pyCMBS provides an easy interface to enable to run any kind of data analysis script. The following minimal requirements need to be fulfilled by the analysis script:

    Programming language: there are no restriction to a programming language. The only requirement is, that the script is stored in an ASCII file.
        Shell: pyCMBS will envoke a shell in which the programm will run. Your analysis programm needs to be able therefore to run on a shell. Some examples are
                Matlab: matlab -nosplash -nodesktop -r <SCRIPTNAME>
                        External python programm: xxxxxxxxxxxxxxxxxxx
                                NCL: xxxxxxxxxxxxxxxxxxxxxxxxx
                                        shell script: no additional information needed
                                            Output directory: The output script needs to write all output to a directory. This directory is supposed to be specified in the programm by a variable.
                                                    The output script is expected to check if its output directory is already existing. If this is not the case, it is expected to generate the output directory by its own. pyCMBS will not know the output directory and can therefore not check for it!
                                                        Graphical output: it is recommended that all graphical output is written to graphic files. Graphical output will in general work also with pyCMBS, but the programm is likely to perform faster and also more robust if every graphical output is written to file.
                                                            Success report: pyCMBS needs to know, if the programm has terminated correctly. The exitcode of the shell will be checked for that purpose. If this is zero, then pyCMBS assumes that everything went fine. However, if the programm has no appropriate error handling, then the return code might be zero even if the programm did not terminate successfully. pyCMBS supports therefore a second way to document successfully programm execution. After completion of the programm, pyCMBS will check for a file with the name <SCRIPTNAME>.log. This file is expected to be a simple ASCII file with a single line. In case of successfull programm completion, the file needs to contain the string "TRUE". Any other entry will be considered as a not successful completion of the programme. It is suggested that your script creates the log file in the beginning of its exectution with a value FALSE, and then overwrites this value with TRUE in the very last line of the code.




xxx::

    program test
    implicit none
    integer F,i

    open(10,file='test.x.log');
    write(10,*) 'FALSE'; close(10)

    do i=1,100
        print*, 'Doing some numbercrunching'
    end do

    !finished successfully here
    open(10,file='test.x.log')
    write(10,*) 'TRUE'; close(10)

    end program

blabla






