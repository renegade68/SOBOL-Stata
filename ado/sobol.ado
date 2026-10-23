*! sobol_mata 0.4.2 - 23oct2026
program define sobol, rclass
    version 16.0
    syntax newvarlist(min=1 max=1024), [N(real 0) START(real 0) SCRamble SEED(string)]
    if ("`scramble'" != "" & "`seed'" == "") | ("`scramble'" == "" & "`seed'" != "") {
        display as error "specify scramble and seed(#) together"
        exit 198
    }
    local mode "none"
    local prng "none"
    local seedtext "."
    if "`scramble'" != "" {
        capture confirm number `seed'
        if _rc {
            display as error "seed() must be an integer in [0,4294967295]"
            exit 198
        }
        if missing(`seed') | `seed'<0 | `seed'>4294967295 | `seed'!=floor(`seed') {
            display as error "seed() must be an integer in [0,4294967295]"
            exit 198
        }
        local seedtext : display %21.0f `seed'
        local seedtext = strtrim("`seedtext'")
        local mode "LMS+shift"
        local prng "MT19937-init_genrand-2002"
    }
    if `n' == 0 local n = _N
    if missing(`n') | `n' < 1 | `n' != floor(`n') | `n' > 2147483647 {
        display as error "n() must be a positive integer no larger than 2147483647"
        exit 198
    }
    if missing(`start') | `start' < 0 | `start' != floor(`start') | `start' + `n' > 2^32 {
        display as error "start() must be a nonnegative integer with start() + n() <= 2^32"
        exit 198
    }
    if _N > 0 & _N != `n' {
        display as error "n() must equal the current number of observations"
        display as error "use an empty dataset to generate a different number of observations"
        exit 198
    }
    local d : word count `varlist'
    sobol_setup
    preserve
    capture noisily {
        if _N == 0 quietly set obs `n'
        foreach v of local varlist {
            quietly generate double `v' = .
        }
        if "`scramble'" == "" {
            mata: _sobol_store(tokens(st_local("varlist")), strtoreal(st_local("start")))
        }
        else {
            mata: _sobol_store(tokens(st_local("varlist")), strtoreal(st_local("start")), strtoreal(st_local("seedtext")))
        }
        local j 0
        local starttext : display %21.0f `start'
        local starttext = strtrim("`starttext'")
        foreach v of local varlist {
            local ++j
            label variable `v' "Sobol coordinate `j' (`mode', Joe-Kuo D6)"
            char `v'[sobol_dimension] "`j'"
            char `v'[sobol_start] "`starttext'"
            char `v'[sobol_bits] "32"
            char `v'[sobol_version] "0.4.2"
            char `v'[sobol_scramble] "`mode'"
            char `v'[sobol_seed] "`seedtext'"
            char `v'[sobol_prng] "`prng'"
        }
    }
    local rc = _rc
    if `rc' {
        restore
        exit `rc'
    }
    restore, not
    if 2^round(ln(`n')/ln(2)) != `n' {
        display as text "note: use n() = 2^m for a complete balanced Sobol net"
    }
    if `start' > 0 {
        display as text "note: start() is for continuation/indexing; arbitrary skipping can impair balance"
    }
    return scalar N = `n'
    return scalar dimension = `d'
    return scalar start = `start'
    return scalar next = `start' + `n'
    return scalar bits = 32
    return scalar seed = `seedtext'
    return local varlist "`varlist'"
    return local version "0.4.2"
    return local scramble "`mode'"
    return local prng "`prng'"
    return local directions "Joe-Kuo D(6), 2010-09-16"
    return local order "Gray code; zero-based index; first point is digital shift (zero if unscrambled)"
end
