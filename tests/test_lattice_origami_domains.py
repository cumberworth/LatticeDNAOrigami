#!/usr/env python

"""Unit tests for the lattice_origami_domains module.

Run tests with py.test.
"""

import json
import sys
import pdb
import pytest
import scipy.constants
from lattice_origami_domains.lattice_origami_domains import *


def all_possible_results_returned_set(result_set, function, *args, max_i=100):
    """Rerun function until all possible results returned or max_i reached.

    For testing stoichastic functions with tractable range. Will fail if max_i
    reached or output outside of range given.
    """
    result_set = set(result_set)
    output_set = set()
    i = 0
    while (result_set != output_set) and (i < max_i):
        output = function(*args)
        assert (output in result_set)
        output_set.add(output)
        i += 1

    assert i < max_i

def all_possible_results_returned_list(result_set, function, sort_key, *args, max_i=10000):
    """Rerun function until all possible results returned or max_i reached.

    For testing stoichastic functions with tractable range. Will fail if max_i
    reached or output outside of range given.
    """
    output_set = []
    i = 0
    while (result_set != output_set) and (i < max_i):
        output = function(*args)
        assert (output in result_set)
        if output not in output_set:
            output_set.append(output)
        i += 1
        result_set.sort(key=sort_key)
        output_set.sort(key=sort_key)

    assert i < max_i


def correct_number_possible_results_returned(num_results, function, *args, max_i=1000000):
    """Rerun function many times and check if number unique results expected."""

    output_set = []
    for i in range(max_i):
        output = function(*args)
        if output not in output_set:
            output_set.append(output)

    pdb.set_trace()
    assert len(output_set) == num_results


@pytest.mark.parametrize('value, multiple, expected', [
    (10, 5, True),
    (11, 5, False),
    (0, 5, True),
    (10, 0, False)])
def test_value_is_multiple(value, multiple, expected):
    test_is_multiple = value_is_multiple(value, multiple)
    assert test_is_multiple == expected


@pytest.mark.parametrize('vector, rotation_axis, direction, expected', [
    ([1, 0, 0], XHAT, 1, [1, 0, 0]),
    (np.array([1, 0, 0]), XHAT, 1, np.array([1, 0, 0])),
    (np.array([3, 7, 2]), XHAT, 1, np.array([3, -2, 7])),
    (np.array([3, 7, 2]), XHAT, -1, np.array([3, 2, -7])),
    (np.array([3, 7, 2]), YHAT, 1, np.array([2, 7, -3])),
    (np.array([3, 7, 2]), ZHAT, 1, np.array([-7, 3, 2]))])
def test_rotate_vector_quarter(vector, rotation_axis, direction, expected):
    test_vector = rotate_vector_quarter(vector, rotation_axis, direction)
    try:
        assert all(test_vector == expected)
    except TypeError:
        assert test_vector == expected


@pytest.mark.parametrize('sequence, T, expected', [

    # Two terminal AT
    ('TTATAACT', 300, -1424.1117624983883),
# (0.2 - 300 * -0.0057 + -7.6 - 300 * -0.0213 + -7.2 - 300 * -0.0204 + -7.6 - 300 * -0.0213 + -7.8 - 300 * -0.021 + 2 * (2.2 - 300 * 0.0069)) * 4.184 * 1000 / scipy.constants.gas_constant

    # Complimentary to previous (and reversed to be 5' to 3')
    ('AGTTATAA', 300, -1424.1117624983883),

    # Palindrome
    ('TGCATGCA', 320, -1252.0106237088303),
# (0.2 - 320 * -0.0057 + 4 * (-8.5 - 320 * -0.0227) + 2 * (2.2 - 320 * 0.0069) - 320 * -0.0014) * 4.184 * 1000 / scipy.constants.gas_constant

    # One terminal AT
    ('CTATAACT', 300, -1635.4640382048626),
# (0.2 - 300 * -0.0057 + -7.8 - 300 * -0.021 + -7.2 - 300 * -0.0204 + -7.6 - 300 * -0.0213 + -7.8 - 300 * -0.021 + 2.2 - 300 * 0.0069) * 4.184 * 1000 / scipy.constants.gas_constant

    # No terminal AT
    ('CTATAACG', 300, -2173.909121552311),
# (0.2 - 300 * -0.0057 + -7.8 - 300 * -0.021 + -7.2 - 300 * -0.0204 + -7.6 - 300 * -0.0213 + -10.6 - 300 * -0.0272) * 4.184 * 1000 / scipy.constants.gas_constant

    # Sequence with every pair in the table
    ('AAATTACAGTCTGACGGCGG', 300, -7256.4281325889615)])
# (0.2 - 300 * -0.0057 -7.6 - 300 * -0.0213 + -7.2 - 300 * -0.0204 + -7.2 - 300 * -0.0213 + -8.5 - 300 * -0.0227 + -8.4 - 300 * -0.0224 + -7.8 - 300 * -0.021 + -8.2 - 300 * -0.0222 + -10.6 - 300 * -0.0272 + -9.8 - 300 * -0.0244 + -8 - 300 * -0.0199 + 2.2 - 300 * 0.0069) * 4.184 * 1000 / scipy.constants.gas_constant
def test_calc_hybridization_energy(sequence, T, expected):
    cation_M = 1
    energy = calc_hybridization_energy(sequence, T, cation_M)
    assert math.isclose(energy, expected)


# calc_complimentary_sequence was implicitly tested with test_calc_hybridization


@pytest.mark.parametrize('sequence, expected', [

    # Regular palindromes should fail
    ('TGAAGT', False),

    # Palindromic sequences are those which, when complimented and reversed, are
    # equal (i.e., they fold back on themselves and hybridize (hairpin))
    ('TGATCA', True)])
def test_sequence_is_palidromic(sequence, expected):
    assert sequence_is_palindromic(sequence) == expected


@pytest.fixture
def example_origami_json():
    return JSONInputFile('example_origami.json')


class TestJSONInputFile:

    def test_identities(self, example_origami_json):
        expected = [[-1, -2, -3, -4, -5, -6], [1, 2], [3, 4], [5, 6]]
        assert example_origami_json.identities == expected

    def test_sequences(self, example_origami_json):
        expected = ['TCCCTAGA', 'GGGTGGGA', 'CTCAAAGG', 'TTGTTGAA', 'GGAATAAG', 'GCTAGCGG']
        assert expected == example_origami_json.sequences

    def test_chains(self, example_origami_json):
        expected = [
            {
                'index': 0,
                'identity': 0,
                'positions': [[0, 0, 0], [1, 0, 0], [2, 0, 0], [3, 0, 0], [4, 0, 0], [5, 0, 0]],
                'orientations': [[0, 1, 0], [0, 0, 1], [0, 1, 0], [0, 0, -1], [0, 0, -1], [0, 0, 1]]
            }, {
                'index': 1,
                'identity': 1,
                'positions': [[1, 0, -1], [1, 0, 0]],
                'orientations': [[0, -1, 0], [0, 0, -1]]
            }, {
                'index': 2,
                'identity': 2,
                'positions': [[2, 0, 0], [3, 0, 0]],
                'orientations': [[0, -1, 0], [0, 0, 1]]
            }
        ]
        chains = example_origami_json.chains(0)
        assert chains == expected


@pytest.fixture
def example_origami_system(example_origami_json):
    return OrigamiSystemEight(example_origami_json, 0, 300, 1)


class TestOrigamiSystemEight:

    def test_chains(self, example_origami_system, example_origami_json):
        assert example_origami_system.chains == example_origami_json.chains(0)

    @pytest.mark.parametrize('chain_index, domain_index, expected', [
        (0, 0, np.array([0, 0, 0])),
        (2, 1, np.array([3, 0, 0]))])
    def test_get_domain_position(self, chain_index, domain_index, expected,
            example_origami_system):
        position = example_origami_system.get_domain_position(chain_index,
                domain_index)
        assert all(position == expected)

    @pytest.mark.parametrize('chain_index, domain_index, expected', [
        (0, 0, np.array([0, 1, 0])),
        (2, 1, np.array([0, 0, 1]))])
    def test_get_domain_orientation(self, chain_index, domain_index, expected,
            example_origami_system):
        orientation = example_origami_system.get_domain_orientation(chain_index,
                domain_index)
        assert all(orientation == expected)

    @pytest.mark.parametrize('position, expected', [
        (np.array([0, 0, 0]), UNBOUND),
        (np.array([1, 0, 0]), BOUND),
        (np.array([4, 0, 0]), UNBOUND)])
    def test_get_position_occupancy(self, position, expected,
            example_origami_system):
        test_occupancy = example_origami_system.get_position_occupancy(position)
        assert test_occupancy == expected

    @pytest.mark.parametrize('chain_index, domain_index, expected', [
        (0, 0, UNBOUND),
        (1, 0, UNBOUND),
        (1, 1, BOUND),
        (0, 1, BOUND)])
    def test_get_domain_occupancy(self, chain_index, domain_index, expected,
            example_origami_system):
        test_occupancy = example_origami_system.get_domain_occupancy(chain_index,
                domain_index)
        assert test_occupancy == expected

    @pytest.mark.parametrize('domain, expected', [
        ((0, 1), (1, 1)),
        ((1, 1), (0, 1)),
        ((0, 4), ())])
    def test_get_bound_domain(self, domain, expected, example_origami_system):
        test_domain = example_origami_system.get_bound_domain(*domain)
        assert test_domain == expected

    def test_get_random_staple_identity(monkeypatch, example_origami_system):
        identities = example_origami_system.identities
        identity_set = set(range(1, len(identities)))
        test_identity_set = set()
        while test_identity_set != identity_set:
            test_identity = example_origami_system.get_random_staple_identity()[0]
            assert (test_identity in identity_set)
            test_identity_set.add(test_identity)

    @pytest.mark.parametrize('domain, sequence', [

        ((0, 0), 'TCCCTAGA'),
        ((1, 1), 'GGGTGGGA')])
    def test_get_hybridization_energy(self, domain, sequence,
                example_origami_system):
        cation_M = 1
        expected = calc_hybridization_energy(sequence, 300, cation_M)
        test_energy = example_origami_system.get_hybridization_energy(*domain)
        assert test_energy == expected

    def test_set_unbound_to_unassigned(self,
            example_origami_system):

        # Move an unbound domain to an unassigned position
        chain_index = 0
        domain_index = 0
        position = np.array([1, 0, 1])
        orientation = np.array([0, 1, 0])
        example_origami_system.unassign_domain(chain_index, domain_index)
        example_origami_system.set_domain_configuration(chain_index,
                domain_index, position, orientation)

    def test_set_matching_unbound_to_unbound_same_helix(self,
            example_origami_system):

        # Move an unbound domain to a bound domain part of a another helix
        # First set correct orientation on static domain
        chain_index = 0
        domain_index = 0
        orientation = np.array([0, -1, 0])
        example_origami_system.set_domain_orientation(chain_index,
                domain_index, orientation)

        chain_index = 1
        domain_index = 0
        position = np.array([0, 0, 0])
        orientation = np.array([0, 1, 0])
        example_origami_system.unassign_domain(chain_index, domain_index)
        example_origami_system.set_domain_configuration(chain_index,
                domain_index, position, orientation)

    def test_set_matching_unbound_to_unbound_new_helix(self,
            example_origami_system):

        # Move an unbound domain to a bound domain contiguous to a helical
        # domain but part of a new helix and obeying twist constraints
        # First set correct staple domain position and orientation
        chain_index = 1
        domain_index = 0
        position = np.array([1, 0, 1])
        orientation = np.array([0, 0, 1])
        example_origami_system.unassign_domain(chain_index, domain_index)
        example_origami_system.set_domain_configuration(chain_index,
                domain_index, position, orientation)

        # Unbind second domain on the staple
        domain = (1, 1)
        position = np.array([1, -1, 1])
        orientation = np.array([1, 0, 0])
        example_origami_system.unassign_domain(*domain)
        example_origami_system.set_domain_configuration(*domain, position,
                orientation)

        # Move scaffold domain
        chain_index = 0
        domain = 0
        position = np.array([1, 0, 1])
        orientation = np.array([0, 0, -1])
        example_origami_system.unassign_domain(chain_index, domain_index)
        example_origami_system.set_domain_configuration(chain_index,
                domain_index, position, orientation)

    def test_bind_non_complimentary_correct_orientation(self,
            example_origami_system):

        # Attempt to bind non-complimentary domains in correct orientation
        # First set correct orientation on static domain
        chain_index = 0
        domain_index = 0
        orientation = np.array([0, -1, 0])
        example_origami_system.set_domain_orientation(chain_index,
                domain_index, orientation)

        # Manually change scaffold domain identity
        example_origami_system.identities[0][0] = 2

        # Bind
        chain_index = 1
        domain_index = 0
        position = np.array([0, 0, 0])
        orientation = np.array([0, 1, 0])
        example_origami_system.unassign_domain(chain_index, domain_index)
        with pytest.raises(ConstraintViolation):
            example_origami_system.set_domain_configuration(chain_index,
                domain_index, position, orientation)

    def test_bind_complimentary_wrong_orientation(self, example_origami_system):

        # Attempt to bind complimentary domains in wrong orientations
        # First set correct orientation on scaffold domain
        chain_index = 0
        domain_index = 0
        orientation = np.array([0, -1, 0])
        example_origami_system.set_domain_orientation(chain_index,
                domain_index, orientation)

        # Iterate through all incorrect staple orientations before binding
        p_system = copy.deepcopy(example_origami_system)
        domain = (1, 0)
        position = np.array([0, 0, 0])
        orientations = [np.array([1, 0, 0]),
                        np.array([-1, 0, 0]),
                       #np.array([0, 1, 0]), correct
                        np.array([0, -1, 0]),
                        np.array([0, 0, 1]),
                        np.array([0, 0, -1])]
        for orientation in orientations:
            p_system.unassign_domain(*domain)
            with pytest.raises(ConstraintViolation):
                p_system.set_domain_configuration(*domain,
                    position, orientation)
            p_system = copy.deepcopy(example_origami_system)

    # Itereate through all incorrect twist orientations before binding (but
    # with correct binding orientation for staple)
        p_system = copy.deepcopy(example_origami_system)
        scaffold_domain = (0, 0)
        staple_domain = (1, 0)
        position = np.array([0, 0, 0])
        orientations = [np.array([1, 0, 0]),
                        np.array([-1, 0, 0]),
                        np.array([0, 1, 0]),
                        #np.array([0, -1, 0]), correct
                        np.array([0, 0, 1]),
                        np.array([0, 0, -1])]
        for orientation in orientations:
            p_system.unassign_domain(*staple_domain)
            p_system.set_domain_orientation(*scaffold_domain,
                    orientation)
            with pytest.raises(ConstraintViolation):
                p_system.set_domain_configuration(*staple_domain,
                    position, -orientation)
            p_system = copy.deepcopy(example_origami_system)

    def test_set_unbound_to_bound(self, example_origami_system):

        # Attempt to set unbound domain to bound position
        domain = (1, 0)
        position = np.array([2, 0, 0])
        orientation = np.array([1, 0, 0])
        with pytest.raises(ConstraintViolation):
            example_origami_system.set_domain_configuration(*domain, position,
                    orientation)

    def test_new_helix_bad_twist(self, example_origami_system):

        # Attempt to set bind domain in new helix contigous to old helix with
        # wrong twist but correct binding orientation and identities.
        # First look at correct position but wrong orientations
        scaffold_domain = (0, 0)
        staple_domain = (1, 0)

        # Staple domain position and orientation
        position = np.array([1, 0, 1])
        orientation = np.array([0, 0, 1])
        example_origami_system.unassign_domain(*staple_domain)
        example_origami_system.set_domain_configuration(*staple_domain, position, orientation)
        p_system = copy.deepcopy(example_origami_system)

        orientations = [np.array([1, 0, 0]),
                        np.array([-1, 0, 0]),
                        np.array([0, 1, 0]),
                        np.array([0, -1, 0])]
                        #np.array([0, 0, 1]) correct
                        #np.array([0, 0, -1])] correct
        for orientation in orientations:
            p_system.set_domain_orientation(*staple_domain,
                    -orientation)
            p_system.unassign_domain(*scaffold_domain)
            with pytest.raises(ConstraintViolation):
                p_system.set_domain_configuration(*scaffold_domain,
                    position, orientation)
            p_system = copy.deepcopy(example_origami_system)

        # Iterate through incorrect positions
#        positions = [np.array([1, 1, 0]), np.array([1, -1, 0]), np.array([1, 0, -1])]
#        for position in positions:
#            example_origami_system.unassign_domain(*staple_domain)
#            example_origami_system.set_domain_configuration(*staple_domain, position, orientation)
#            p_system = copy.deepcopy(example_origami_system)
#            orientations = [np.array([1, 0, 0]),
#                            np.array([-1, 0, 0]),
#                            np.array([0, 1, 0]),
#                            np.array([0, -1, 0]),
#                            np.array([0, 0, 1]),
#                            np.array([0, 0, -1])]
#            for orientation in orientations:
#                p_system.set_domain_orientation(*staple_domain,
#                        -orientation)
#                p_system.unassign_domain(*scaffold_domain)
#                with pytest.raises(ConstraintViolation):
#                    p_system.set_domain_configuration(*scaffold_domain,
#                        position, orientation)
#                p_system = copy.deepcopy(example_origami_system)

    # Consider repeating tests with staple in double bound state, or with one
    # not contiguous to another bound domain domain

    @pytest.mark.parametrize('chain_index, domain_index, expected', [

        # Change unbound scaffold strand domain
        (0, 0, np.array([1, 0, 0])),

        # Change staple strand domain
        (1, 0, np.array([0, 1, 0]))])
    def test_set_domain_orientation_correct(self, chain_index, domain_index, expected,
            example_origami_system):
        example_origami_system.set_domain_orientation(chain_index, domain_index,
                expected)
        test_orientation = example_origami_system.get_domain_orientation(chain_index,
                domain_index)
        assert all(test_orientation == expected)

    @pytest.mark.parametrize('chain_index, domain_index, expected', [

        # Change bound scaffold strand domain
        (0, 1, np.array([1, 0, 0])),

        # Multiply by scalar
        (0, 1, np.array([0, 0, 99]))])
    def test_set_domain_orientation_wrong(self, chain_index, domain_index, expected,
            example_origami_system):
        with pytest.raises(ConstraintViolation):
            example_origami_system.set_domain_orientation(chain_index,
                    domain_index, expected)

    def test_unassign_bound_domain(self, example_origami_system):
        domain = (0, 1)
        position = (1, 0, 0)
        bound_domain = example_origami_system.get_bound_domain(*domain)
        ee_delta = -example_origami_system.get_hybridization_energy(*domain)
        e_delta = example_origami_system.unassign_domain(*domain)
        assert e_delta == ee_delta
        ep_occupancy = UNBOUND
        p_occupancy = example_origami_system.get_position_occupancy(position)
        assert p_occupancy == ep_occupancy
        ed_occupancy = UNASSIGNED
        d_occupancy = example_origami_system.get_domain_occupancy(*domain)
        assert d_occupancy == ed_occupancy
        eb_domain = ()
        b_domain = example_origami_system.get_bound_domain(*domain)
        assert b_domain == eb_domain
        ebd_occupancy = UNBOUND
        bd_occupany = example_origami_system.get_domain_occupancy(*bound_domain)
        assert bd_occupany == ebd_occupancy
        assert example_origami_system._unbound_domains[position] == bound_domain
        
    def test_unassign_unbound_domain(self, example_origami_system):
        domain = (0, 0)
        position = (0, 0, 0)
        example_origami_system.unassign_domain(*domain)
        ep_occupancy = UNASSIGNED
        p_occupancy = example_origami_system.get_position_occupancy(position)
        assert p_occupancy == ep_occupancy
        ed_occupancy = UNASSIGNED
        d_occupancy = example_origami_system.get_domain_occupancy(*domain)
        assert d_occupancy == ed_occupancy
        with pytest.raises(KeyError):
            example_origami_system._unbound_domains[position]

    def test_add_chain(self, example_origami_system):
        identity = 2
        example_origami_system.add_chain(identity)
        expected_unique_indices = [0, 1, 2, 3]
        assert example_origami_system._working_to_unique == expected_unique_indices
        test_position = example_origami_system.get_domain_position(3, 0)
        assert test_position == []
        test_occupancy = example_origami_system.get_domain_position(3, 0)
        assert test_occupancy == []

    def test_delete_chain(self, example_origami_system):
        chain_index = 1

        position_unbound = (1, 0, -1)
        position_bound = (1, 0, 0)

        ee_delta = -example_origami_system.get_hybridization_energy(1, 1)
        e_delta = example_origami_system.delete_chain(1)
        assert e_delta == ee_delta

        e_num_chains = 2
        assert len(example_origami_system._working_to_unique) == e_num_chains
        assert len(example_origami_system._chain_identities) == e_num_chains
        assert len(example_origami_system._positions) == e_num_chains
        assert len(example_origami_system._orientations) == e_num_chains
        assert len(example_origami_system.chain_lengths) == e_num_chains

        p_occupancy = example_origami_system.get_position_occupancy(position_bound)
        assert p_occupancy == UNBOUND
        eup_occupancy = example_origami_system.get_position_occupancy(position_unbound)
        assert eup_occupancy == UNASSIGNED
        bd_occupancy = example_origami_system.get_domain_occupancy(0, 1)
        assert bd_occupancy == UNBOUND
        assert example_origami_system._unbound_domains[position_bound] == (0, 1)

        assert example_origami_system.get_domain_occupancy(1, 0) == BOUND
        assert example_origami_system.get_domain_occupancy(1, 1) == BOUND

        ne_delta = example_origami_system.delete_chain(1)
        assert ne_delta != e_delta

    def test_center(self, example_origami_system):
        example_origami_system.center()

    def test_domains_match(self, example_origami_system):

        # Complimentary domains in correct orientation
        domain_1 = (0, 1)
        domain_2 = (1, 1)
        test_match = example_origami_system._domains_match(*domain_1, *domain_2)
        assert test_match == True

        # Complimentary domains in incorrect orientation
        domain_2_orientation = example_origami_system.get_domain_orientation(
                *domain_2)
        domain_2_orientation = rotate_vector_quarter(domain_2_orientation, XHAT, 1)
        example_origami_system._orientations[domain_2[0]][domain_2[1]] = domain_2_orientation
        test_match = example_origami_system._domains_match(*domain_1, *domain_2)
        assert test_match == False

        # Non-complimentary domains in correct orientation

        # First double check revert orientation goes back to a match
        domain_2_orientation = rotate_vector_quarter(domain_2_orientation, XHAT, -1)
        example_origami_system._orientations[domain_2[0]][domain_2[1]] = domain_2_orientation
        test_match = example_origami_system._domains_match(*domain_1, *domain_2)
        assert test_match == True

        # Delete chain and replace with new chain with different identity
#        example_origami_system.delete_chain(1)
#        example_origami_system.add_chain(3)
#        domain_2 = (2, 0)
#        test_match = example_origami_system._domains_match(*domain_1, *domain_2)
#        assert test_match == False

        # Double check that directly reverting chain identity will return match
#        example_origami_system._chain_identities[2] = 1
#        test_match = example_origami_system._domains_match(*domain_1, *domain_2)
#        assert test_match == True


class TestIO:

    def test_read_write_json_configuration(self, example_origami_system):
        filename = 'json_test.json'
        output_file = JSONOutputFile(filename, example_origami_system, 1)
        output_file._write_configuration(example_origami_system, 0)
        test_input = JSONInputFile(filename)
        assert test_input.identities == example_origami_system.identities
        assert test_input.sequences == example_origami_system.sequences
        assert test_input.chains(0) == example_origami_system.chains

    def test_read_write_hdf5_configuration(self, example_origami_system):
        filename = 'hdf5_test.hdf5'
        output_file = HDF5OutputFile(filename, example_origami_system, 1)
        output_file._write_configuration(example_origami_system, 0)
        test_input = HDF5InputFile(filename)
        assert test_input.identities == example_origami_system.identities
        assert test_input.sequences == example_origami_system.sequences
        assert test_input.chains(0) == example_origami_system.chains

    # Should add tests for adding and removing chains at different timesteps


@pytest.fixture
def test_sim(example_origami_system, example_origami_json):
    move_settings = {MOVETYPE.INSERT_STAPLE: 0.2,
                     MOVETYPE.DELETE_STAPLE: 0.2,
                     MOVETYPE.REGROW_STAPLE: 0.2,
                     MOVETYPE.REGROW_SCAFFOLD_AND_BOUND_STAPLES: 0.2,
                     MOVETYPE.ROTATE_ORIENTATION_VECTOR: 0.2}
    return GCMCBoundStaplesSimulation(example_origami_system,
            move_settings, example_origami_json)


@pytest.fixture
def simple_loop_json():
    return JSONInputFile('simple_loop.json')


@pytest.fixture
def simple_loop(simple_loop_json):
    return OrigamiSystemEight(simple_loop_json, 0, 300, 1)


@pytest.fixture
def test_sim_loop(simple_loop, simple_loop_json):
    move_settings = {MOVETYPE.INSERT_STAPLE: 0.2,
                     MOVETYPE.DELETE_STAPLE: 0.2,
                     MOVETYPE.REGROW_STAPLE: 0.2,
                     MOVETYPE.REGROW_SCAFFOLD_AND_BOUND_STAPLES: 0.2,
                     MOVETYPE.ROTATE_ORIENTATION_VECTOR: 0.2}
    return GCMCBoundStaplesSimulation(simple_loop,
            move_settings, simple_loop_json)


class TestGCMCBoundStaplesSimulation:

    def test_select_movetype(self, test_sim):

        movetypes = set([test_sim._insert_staple, test_sim._delete_staple,
                test_sim._regrow_staple, test_sim._regrow_scaffold_and_bound_staples,
                test_sim._rotate_orientation_vector])
        all_possible_results_returned_set(movetypes, test_sim._select_movetype)

    @pytest.mark.parametrize('ratio, rand_ratio, expected', [
        (2, 0, True),
        (0.5, 0.5, False),
        (0.5, 0.4, True),
        (0.5, 0.6, False)])
    def test_test_acceptance(self, ratio, rand_ratio, expected, test_sim, monkeypatch):
        monkeypatch.setattr(random, 'random', lambda: rand_ratio)
        assert test_sim._test_acceptance(ratio) == expected

    @pytest.mark.parametrize('delta_e, rand_ratio, expected', [
        (0, 0, True),
        (300, 0.5, False),
        (-300, 0, True)])
    def test_configuration_accepted(self, delta_e, rand_ratio, expected, test_sim, monkeypatch):
        monkeypatch.setattr(random, 'random', lambda: rand_ratio)
        test_sim._delta_e = delta_e
        assert test_sim._configuration_accepted() == expected

    # Maybe test staple insertion/deletion when sure about reasonable number densities

    def test_insert_staple_wrong_binding(self, test_sim, monkeypatch):

        # Attempt to grow staple into constraint violation
        original_system = copy.deepcopy(test_sim._accepted_system)
        staple_domain = 0
        scaffold_domain = 4
        p_dimension = XHAT
        p_direction = 1
        o_dimension = XHAT
        o_direction = 1
        randranges = (staple_domain, scaffold_domain, p_direction, o_direction)
        choices = (p_dimension, o_dimension)
        randrange_output = (i for i in randranges)
        choice_output = (i for i in choices)
        test_sim._accepted_system.get_random_staple_identity = lambda: (2, 0)
        test_sim._staple_insertion_accepted = lambda: True
        monkeypatch.setattr(random, 'randrange', lambda *args: next(randrange_output))
        monkeypatch.setattr(random, 'choice', lambda *args: next(choice_output))
        test_sim._insert_staple()
        assert test_sim._accepted_system.chains == original_system.chains

    def test_insert_staple_correct_double_binding(self, test_sim, monkeypatch):

        # Attempt to grow staple to singly bound state
        # Setup initial and expected accepted system
        expected_system = copy.deepcopy(test_sim._accepted_system)
        expected_system.add_chain(3)
        position = np.array([4, 0, 0])
        orientation = np.array([0, 1, 0])
        expected_system.set_domain_orientation(0, 4, -orientation)
        test_sim._trial_system.set_domain_orientation(0, 4, -orientation)
        expected_system.set_domain_configuration(3, 0, position, orientation)
        position = np.array([5, 0, 0])
        orientation = np.array([0, 0, -1])
        expected_system.set_domain_orientation(0, 5, -orientation)
        test_sim._trial_system.set_domain_orientation(0, 5, -orientation)
        expected_system.set_domain_configuration(3, 1, position, orientation)

        staple_domain = 0
        scaffold_domain = 4
        p_dimension = XHAT
        p_direction = 1
        o_dimension = ZHAT
        o_direction = -1
        randranges = (staple_domain, scaffold_domain, p_direction, o_direction)
        choices = (p_dimension, o_dimension)
        randrange_output = (i for i in randranges)
        choice_output = (i for i in choices)
        test_sim._accepted_system.get_random_staple_identity = lambda: (3, 0)
        test_sim._staple_insertion_accepted = lambda: True
        monkeypatch.setattr(random, 'randrange', lambda *args: next(randrange_output))
        monkeypatch.setattr(random, 'choice', lambda *args: next(choice_output))
        test_sim._insert_staple()
        assert test_sim._accepted_system.chains == expected_system.chains

    def test_delete_staple(self, test_sim):
        all_staples = test_sim._accepted_system
        deleted_staple_1 = copy.deepcopy(test_sim._accepted_system)
        deleted_staple_1.delete_chain(1)
        deleted_staple_2 = copy.deepcopy(test_sim._accepted_system)
        deleted_staple_2.delete_chain(2)
        expected_chains = [all_staples.chains, deleted_staple_1.chains, deleted_staple_2.chains]
        def test_delete_staple(test_sim):
            test_sim = copy.deepcopy(test_sim)
            test_sim._delete_staple()
            return test_sim._accepted_system.chains

        def sort_key(chains):
            return chains[1]['identity']

        all_possible_results_returned_list(expected_chains, test_delete_staple, sort_key, test_sim)

    def test_regrow_staple(self, test_sim_loop, monkeypatch):

        # Going to check if all possible regrowths occur for staple one
        # Create list of arrays containing the position and orientation lists
        possible_configs = [([[0, 0, 0], [0, 0, 1]], [[0, 0, 1], [0, 1, 0]])]

        position_1 = np.array([0, 0, 0])
        for position_2 in [[-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, -1]]:
            for dim in range(3):
                for dir in [1, -1]:
                    orientations = np.array([[0, 0, 1], [0, 0, 0]])
                    orientations[1][dim] = dir
                    positions = np.array([position_1, position_2])
                    possible_configs.append((positions.tolist(), orientations.tolist()))

        position_2 = np.array([0, 0, 1])
        for position_1 in [[-1, 0, 1], [0, 1, 1], [0, -1, 1], [0, 0, 2]]:
            for dim in range(3):
                for dir in [1, -1]:
                    orientations = np.array([[0, 0, 0], [0, 1, 0]])
                    orientations[0][dim] = dir
                    positions = np.array([position_1, position_2])
                    possible_configs.append((positions.tolist(), orientations.tolist()))

        def test_regrow_staple(test_sim_loop):
            test_sim = copy.deepcopy(test_sim_loop)
            test_sim._regrow_staple()
            positions = test_sim._accepted_system.chains[1]['positions']
            orientations = test_sim._accepted_system.chains[1]['orientations']
            return (positions, orientations)

        def sort_key(staple_config):
            return staple_config

        #def randrange_gen_maker():
        #    i = 0
        #    while True:
        #        if i == 0:
        #            yield lambda *args: 1
        #            #yield 1
        #        else:
        #            yield random.randrange
#
#                i += 1
#
        #def call_next(*args):
        #    rand_func = next(randrange_poo)
        #    output = rand_func(*args)
        #    return output
#
#        randrange_poo = randrange_gen_maker()
        # For some reason the generator raises StopIteration when I try and do
        # it the way above
        poopoo = [lambda *args: 1] + [random.randrange] * 100000
        test_sim_loop._configuration_accepted = lambda *args: True
        randrange_gen = (i for i in poopoo)
        monkeypatch.setattr(random, 'randrange', lambda *args: next(randrange_gen)(*args))
        all_possible_results_returned_list(possible_configs, test_regrow_staple, sort_key, test_sim_loop)

#    def test_regrow_scaffold_and_bound_staples(self, test_sim_loop, monkeypatch):
#        """This test sucks. It requires that I change the source code. But I
#        don't see how else I can change only specific instances of randrange."""
#
#        # Attempt to regrow first domain and bound staple.
#        #start_domain_index = 0
#        #direction = 0
#        total_configs = 205
#
#        #poopoo = [lambda *args: 0] * 2 + [random.randrange] * 100000000
#        #test_sim_loop._configuration_accepted = lambda *args: True
#        #randrange_gen = (i for i in poopoo)
#        #monkeypatch.setattr(random, 'randrange', lambda *args: next(randrange_gen)(*args))
#
#        def test_regrow_scaffold(test_sim_loop):
#            test_sim = copy.deepcopy(test_sim_loop)
#            test_sim._regrow_scaffold_and_bound_staples()
#            staple_positions = test_sim._accepted_system.chains[1]['positions']
#            staple_orientations = test_sim._accepted_system.chains[1]['orientations']
#            scaffold_position = test_sim._accepted_system.chains[0]['positions'][0]
#            scaffold_orientation = test_sim._accepted_system.chains[0]['orientations'][0]
#            return ((scaffold_position, scaffold_orientation), (staple_positions, staple_orientations))
#
#        correct_number_possible_results_returned(total_configs, test_regrow_scaffold, test_sim_loop)

    def test_rotate_orientation_vector(self, test_sim):
        test_sim._rotate_orientation_vector()
