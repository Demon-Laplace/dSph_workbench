# Current spatial-member video batch

The current entry points are maintained centrally in Code/simulation/server_video_repair and deployed as a pinned independent job package. One mesopsl1/def/mesopsl1_def_long allocation uses one node, 16 CPUs and sequential model processing. This revision covers the 94 available production runs (91 model directories, including orbital branches); it does not change ICs or restart evolution.

Use current_spatial_properties_seeded_v2, astropy_native_direct_v1 and heliocentric_tan_xi_eta_2p1_v1. Source IDs do not select members, prove origin or reject physical rows. Repeated IDs remain rows. Recorded stellar ages/formation times are independent of IDs. Gas uses all local TP0 particles within strict 3D5 kpc, loose3D15 kpc and the5–15 annulus, with matched-aperture HI fractions.

Existing measured dwarf/host trajectories supply only search seeds. Each frame remeasures local centres from current physical stars with the shared spatial helper. Host bulk remains the mass mean inside true3D5 kpc. The standalone package supplies physical_helpers.py and spatial_context.py through DSPH_SPATIAL_HELPERS_DIRECTORY and a slim per-run seed context through DSPH_SPATIAL_CONTEXT.

DSPH_MEDIA_REUSE=scientific explicitly selects ordinary resume by observer/member/display/render-revision identity and existing outputs. It avoids source/file hashes, per-frame certification and post-write rereads; real I/O, physical computation and ffmpeg errors remain errors. Old frozen-ID CSVs/frames are not relabelled as current spatial products. Fully matching caches may resume, and fully matching frames may use --video-only.

Published videos record source/run identity, snapshot numbers, observer, member-selection and display definitions. Old videos and frozen releases remain historical inputs; this batch does not delete them. Submit the single batch and stop: no Slurm monitoring or startup-validation loop is required by the owner. Code has not been exercised on real snapshots in this maintenance step; successful submission is not a completed video claim.
