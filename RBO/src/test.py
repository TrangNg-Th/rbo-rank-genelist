"""
This module contains some test cases.
"""



import string
from rbo_origin import rbo
from rbo_modif import rbo_modified

TESTS = [
    # Sanity checks
    (string.ascii_lowercase, string.ascii_lowercase, 1.0),
    (string.ascii_lowercase, string.ascii_lowercase[:7], 1.0),
    ("abcde", "fghij", 0.0),

    # RBO Paper Figure 5
    ("abcdefg", "zcavwxy", 0.312),

    # Source:  https://ragrawal.wordpress.com/2013/01/18/comparing-ranked-list/
    ("abcde", "bacde", 0.8),
    ("abcde", "abced", 0.95),

    # One-Element lists
    ("a", "a", 1.0),
    ("a", "b", 0),

    # Empty lists
    ("", "", 1),
    ("a", "", 0),
    ("", "a", 0),
]


def test_rbo(list_1: list, list_2: list, expected: float):
    """
    Args:
        list_1: List 1.
        list_2: List 2.
        expected: Expected RBO.
    Returns:
        None
    """
    p = 0.95  # pylint: disable=invalid-name
    list_1, list_2 = list(list_1), list(list_2)
    print("List 1 is: {}".format(list_1))
    print("List 2 is: {}".format(list_2))

    rbo_res = rbo(list_1, list_2)
    print(f"The implemented Average Overlap is: {rbo_res}")
    print(f"The correct answer is:              {expected}")
    assert np.round(rbo_res, decimals=3) == expected, "Not correct"
    
    
# TEST CASE FOR MODIFIED RBO WITH SAME WEIGHT FUNCTION
# Geometrical progression
def wg_geom(d, p=0.5):
    """
    Arguments:
        d -- rank of eval
        p -- parameter of geometrical sum
    return the geometical progression up to rank d with parameter p
    
    """
    return(1.0 * (1 - p) * p**d)
    
    
    
def test_rbo_vs_rbo_modified(list_1: list, list_2: list, expected: float):
    """
    Args:
        list_1: List 1.
        list_2: List 2.
        expected: Expected RBO.
    Returns:
        None
    """
    p = 0.95  # pylint: disable=invalid-name
    list_1, list_2 = list(list_1), list(list_2)
    print("List 1 is: {}".format(list_1))
    print("List 2 is: {}".format(list_2))

    rbo_res = rbo(list_1, list_2, p=0.5)
    rbo_mod = rbo_modified(wg_geom, list_1, list_2)
    print(f"The result of rbo: {rbo_res}")
    print(f"The result of rbo_modified: {rbo_mod}")
 
    assert np.round(rbo_res, decimals=3) == np.round(rbo_mod, decimals=3), \
    "The results are not equal for this test case, please check back rbo_modified"
    
 