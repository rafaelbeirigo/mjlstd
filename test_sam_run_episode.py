import sam_run_episode
import unittest
import numpy as np
import numpy.testing as npt


class TestSamRunEpisode(unittest.TestCase):
    def test_get_next_theta_works_for_one_dimension(self):
        P = np.array([[1.]])
        theta = 0
        changed = False
        for _ in range(1024):
            if sam_run_episode.get_next_theta(theta, P) != 0:
                changed = True
                break

        self.assertFalse(changed)

    def test_get_next_theta_works_for_eye_two(self):
        P = np.eye(2)
        for theta in range(2):
            changed = False
            for _ in range(1024):
                if sam_run_episode.get_next_theta(theta, P) != theta:
                    changed = True
                    break
            if changed:
                break

        self.assertFalse(changed)

    @unittest.skip("Too long to run all the time.")
    def test_get_next_theta_works_for_uniform_three(self):
        N = 3
        P = (1. / N) * np.ones((N, N))
        num_samples = 10240
        P_sample = np.zeros_like(P)
        for theta in range(N):
            for _ in range(num_samples):
                next_theta = sam_run_episode.get_next_theta(theta, P)
                P_sample[theta, next_theta] += 1
                theta = next_theta

        P_sample /= num_samples

        npt.assert_array_almost_equal(P, P_sample, decimal=2)

    def test_get_next_theta_rejects_non_row_stochastic_matrices(self):
        self.assertRaises(ValueError,
                          sam_run_episode.get_next_theta,
                          0,
                          np.zeros((2, 2)))

    def test_get_next_theta_rejects_non_square_matrices(self):
        self.assertRaises(ValueError,
                          sam_run_episode.get_next_theta,
                          0,
                          np.zeros((2, 2)))

    def test_get_next_theta_accepts_only_non_negative_matrices(self):
        self.assertRaises(ValueError,
                          sam_run_episode.get_next_theta,
                          0,
                          -np.ones((2, 2)))
