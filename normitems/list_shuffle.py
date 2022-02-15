import random
import numpy as np

def shuffle_list(some_list):

    randomized_list = some_list[:]
    while True:
        random.shuffle(randomized_list)
        for a, b in zip(some_list, randomized_list):
            if a == b:
                break
        return randomized_list

def find_unique_shuffle_lists(some_list, nshuffle):

    list_of_randomized_lists = []
    ngood = 0

    while True:
        randomized_list = shuffle_list(some_list)
        good = True

        for previous_list in list_of_randomized_lists:

            if previous_list == randomized_list:
                good = False
                break
        if good:
            ngood += 1
            list_of_randomized_lists.append(randomized_list)

        if ngood == nshuffle:
            break

    return list_of_randomized_lists

if __name__ == "__main__":

    a_list = list( np.arange(10) )
    print( a_list )
    print( shuffle_list(a_list) )
    nshuffle = 10
    print( find_unique_shuffle_lists(a_list, nshuffle) )
