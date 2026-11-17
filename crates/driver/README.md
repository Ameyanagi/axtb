# axtb-driver

The request/response layer of the [axtb](https://github.com/atomsforhumanity/axtb) front ends:
the command-line program ([`axtb-cli`](https://crates.io/crates/axtb-cli)) and the Python
package (`axtb` on PyPI). A `Request` names the method (GFN2-xTB, g-xTB v1.1 or 2.0.1,
GFN-FF), the structure and the settings; `run` checks it, dispatches it, and returns a
`Response` with the same fields for every front end, including Hessians and their vibrational
analysis. A stateful `Calculator` runs one geometry after another.

Every response names the numerics identity of the code that produced it, `numerics`, which
`numerics(method, composite)` gives before a run: `"<method>:<epoch>"` (`"gfn2:1"`,
`"gxtb-v1.1:1"`, `"gxtb-v2.0.1:1"`, `"gfnff:1"`; the g-xTB solvation composite
`"gxtb-v2.0.1:1+gfn2:1"`). One epoch per method covers all its conventions; it changes with every
commit that may change that method's results (round-off or model changes, R or M in the
changelog) and never with one that keeps them bit for bit. The history of the epochs is in the
changelog.

License: LGPL-3.0-or-later (`COPYING.LESSER`, with `COPYING`). `NOTICE` credits the upstream
work.
