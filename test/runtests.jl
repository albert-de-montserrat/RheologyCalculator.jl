function runtests()
    files = readdir(@__DIR__)
    test_files = filter(startswith("test_"), files)

    for f in test_files
        !isdir(f) && include(f)
    end
    return
end

runtests()