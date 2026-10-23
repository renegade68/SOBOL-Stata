*! sobol_mata 0.4.2 - 23oct2026
program define sobol_setup
    version 16.0
    syntax
    local loaded 0
    mata: st_local("loaded", strofreal(findexternal("sobol_version()")!=NULL))
    if `loaded' {
        mata: st_local("installed", sobol_version())
        if "`installed'" != "0.4.2" {
            display as error "another Sobol version is loaded; restart Stata before loading this version"
            exit 498
        }
        exit
    }
    local conflict 0
    mata: st_local("conflict", strofreal(findexternal("sobol()")!=NULL))
    if `conflict' {
        display as error "a different Mata sobol() function is already loaded"
        exit 498
    }
    quietly findfile _sobol_direction_data.mata
    local datafile "`r(fn)'"
    quietly findfile sobol.mata
    local sourcefile "`r(fn)'"
    quietly findfile _sobol_mt19937.mata
    local rngfile "`r(fn)'"
    quietly findfile sobol_rqmc.mata
    local rqmcfile "`r(fn)'"
    quietly findfile sobol_sa.mata
    local safile "`r(fn)'"
    quietly run "`datafile'"
    quietly run "`rngfile'"
    quietly run "`sourcefile'"
    quietly run "`rqmcfile'"
    quietly run "`safile'"
end
