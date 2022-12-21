function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--Decomp"
            help = "Domain decomposition method"
            arg_type = String
            default = "Hilbert"
        "--ending"
            help = "ending value of array"
            arg_type = Int
            default = 10
        "--step"
            help = "step value of array"
            arg_type = Int
            default = 1
    end
    return parse_args(s)
end
