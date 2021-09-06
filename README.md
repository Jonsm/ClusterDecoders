#  Readme

Utilities for testing a decoder for the biased color code. This is the 2D color code with a local unitary applied to every other site (i.e. each site in one of the bipartitions of the sublattice). The unitary takes X->Y->iZ. This code has the property that optimal correction gives a 50% threshold at 100% Z bias (only Z errors).

Here, I develop a local decoder based on the RG decoder of https://arxiv.org/pdf/1112.3252.pdf. This decoder can correct at arbitrary bias.

The main decoder class is cluster_decoder. Utilities for testing it are ca_sim, thermal_bias_sim, and threshold_sim.
