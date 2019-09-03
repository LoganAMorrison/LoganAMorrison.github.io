using Weave

parent_dir = string(@__DIR__)
dirs = filter(x -> isdir(joinpath(parent_dir, x)), readdir(parent_dir))

for dir in dirs
    if dir[1] == '_'
        files = readdir(parent_dir * "/" * dir)
        weave(parent_dir * "/" * dir * "/" * files[1], out_path=parent_dir)
    end
end
