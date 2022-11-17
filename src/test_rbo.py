# THIS MODULE CONTAINS TESTS FOR RBO MEASURE

import numpy as np
import string
from rbo import rbo
from rbomod import


TESTS = [
    # Sanity checks
    (string.ascii_lowercase, string.ascii_lowercase),
    (string.ascii_lowercase, string.ascii_lowercase[:7]),
    ("abcde", "fghij"),

    # RBO Paper Figure 5
    ("abcdefg", "zcavwxy"),

    # Source:  https://ragrawal.wordpress.com/2013/01/18/comparing-ranked-list/
    ("abcde", "bacde"),
    ("abcde", "abced"),

    # One-Element lists
    ("a", "a"),
    ("a", "b"),

    # Empty lists
    ("", ""),
    ("a", ""),
    ("", "a"),
]


def test_rbo(list_1: list, list_2: list, p=0.5):
    """
    Args:
        list_1: List 1.
        list_2: List 2.
        expected: Expected RBO.
    Returns:
        None
    """
  
    list_1, list_2 = list(list_1), list(list_2)
    print('-'*50)
    print("List 1 is: {}".format(list_1))
    print("List 2 is: {}".format(list_2))

    rbo_res = rbo(list_1, list_2, p=p)
    print(f'The chosen p is {p}')
    print(f"The implemented Average Overlap is: {rbo_res}")
    
    
    
# TEST CASE FOR MODIFIED RBO WITH SAME WEIGHT FUNCTION     
def test_rbo_vs_rbo_modified(list_1: list, list_2: list, p=0.5):
    """
    Args:
        list_1: List 1.
        list_2: List 2.
        expected: Expected RBO.
    Returns:
        None
    """
   
    list_1, list_2 = list(list_1), list(list_2)
    print('-'*50)
    print("List 1 is: {}".format(list_1))
    print("List 2 is: {}".format(list_2))
    
    rbo_res = rbo(list_1, list_2, p=p)
    rbo_mod = rbo_modified(lstS=list_1, lstT=list_2, wg_func=wg_geom, params_wg_func={'geometric':None, 'p':p})
    print(f'The chosen p is {p}')
    print(f"The implemented Average Overlap is: {rbo_res}")
    print(f"The result of rbo_modified:         {rbo_mod}")
    assert np.round(rbo_res, decimals=3) == np.round(rbo_mod, decimals=3), \
    "The results are not equal for this test case, please check back rbo_modified"
    
   
    
for p in [k*0.1 for k in range(1, 10)]:
    print('*'*80)
    for test in TESTS:
        test_rbo(test[0], test[1], p=p)
        test_rbo_vs_rbo_modified(list_1=test[0], list_2=test[1], p=p)