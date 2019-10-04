using Weave
dir = string(@__DIR__)
weave_dir = joinpath(dir, "weave_source")
out_dir = joinpath(dir, "_posts")
fig_dir = joinpath("..", "assets", "img")
filesnames = ["2019-08-21-thermal-distribution-functions.jmd",
               "2019-08-22-thermal-cross-section-eta-prime-delta.jmd",
               "2019-08-23-thermal-cross-section-2-to-4-eta-prime.jmd",
               "2019-08-24-one-loop-effective-potential-in-the-abelian-higgs-and-gauge-dependence.jmd",
               "2019-08-24-solving-the-boltzmann-equation.jmd",
               "2019-08-25-temperature-evolution-of-a-decoupled-species.jmd"]

for file_name in files_names
    weave(joinpath(weave_dir, file_name),
          out_path=out_dir,
          doctype="github",
          fig_path=fig_dir)
end

run(`python fix_image_location.py $filesnames`)
