# Reviewer-derived corollary: the complete three-copy qutrit threshold

This is a consequence of the submitted endpoint theorem, not a claim that the
original manuscript already stated the corollary. It needs no new SDP solve.

Let rho_alpha = (I + alpha F)/(9 + 3 alpha), -1 <= alpha <= 1. Assuming the
submitted endpoint theorem, rho_alpha is three-copy distillable if and only if
alpha < -1/2.

## Proof

Use the unnormalised maximally entangled vector Omega = sum_{j=0}^2 |jj>, and
write W_alpha = I + alpha |Omega><Omega|. The endpoint theorem says
W_{-1/2}^{tensor 3} is nonnegative on every Schmidt-rank-at-most-two vector.

It also implies the corresponding statement for one and two copies. Indeed,
tensor any lower-copy test vector with |0>_A |1>_B on each extra copy. Each such
factor has expectation one against W_{-1/2}, and Schmidt rank is unchanged.

For -1/2 <= alpha <= 0 put s = -2 alpha, so 0 <= s <= 1 and
W_alpha = s W_{-1/2} + (1-s) I. Expand its third tensor power. Every coefficient
is nonnegative. Each summand consists of k endpoint factors and 3-k identity
factors, where 0 <= k <= 3. Resolve the identity factors in product bases across
Alice and Bob. Every contracted vector still has Schmidt rank at most two, so
every expectation in each summand is nonnegative by the lower-copy statements
just established. For alpha >= 0, W_alpha itself is positive semidefinite.

For alpha < -1/2 use the unnormalised test vector
(|00> + |11>) on the first copy and |01> on each of the other two copies. Its
Schmidt rank across all Alice factors versus all Bob factors is two, and its
expectation against W_alpha^{tensor 3} is 2 + 4 alpha < 0. The normalising
factor (9+3 alpha)^{-3} is positive throughout the state parameter range.

Thus the stated threshold follows. This proof does not establish anything
about four or more copies of the endpoint, or unrestricted higher dimensions.
