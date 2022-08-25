import glob
import os
import sys

### === special string files query
def find_str_func( Path, file_type, target_str, out_file = None):
    """
    Path : where to find files
    --------------------------
    file_type : which type file to be found
    '*.py' : .py files
    '*.o'  : C files
    you can define any type
    --------------------------
    target_str : the special string you want to find out
    like : 'follow_me', 'where_are'... Please note it is 'A' or 'a'
    --------------------------
    out_file : file to list the file name after finding them
    default is None, which means print in the screen or terminal directly
    """

    files = glob.glob( Path + file_type )

    N_file = len( files )

    tag_str = target_str

    out_str = []

    for ll in range( N_file ):

        dat = open( files[ ll ] )

        lines = dat.readlines()
        N_lines = len( lines )

        for dd in range( N_lines ):

            if tag_str in lines[ dd ]:

                out_str.append( files[ ll ] )

                break

            else:
                continue
    if len( out_str ) == 0:
        print('No such files!')
        return

    else:
        if out_file is None:
            print( out_str )
            return

        else:

            doc = open( out_file, 'w')

            nn = len( out_str )

            for dd in range( nn ):

                s = out_str[ dd ]

                print( s, file = doc )

            doc.close()

            return


### === test
if __name__ == "__main__":

    path = '/home/xkchen/mywork/CSST/'
    file_ty = '*.py'
    tagt_str = 'ICL'
    out_file = '/home/xkchen/file_list.txt'

    find_str_func( path, file_ty, tagt_str, out_file = out_file )

