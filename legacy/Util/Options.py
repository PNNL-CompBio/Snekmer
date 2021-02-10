#!/usr/bin/env python

import os, sys, string
from getopt import getopt, GetoptError

"""
A module based on the standard getopt which will handle command line input
based on a specification dictionary provided by the app. The dictionary is
indexed by a descriptive keyword and contains a Tuple with the following
organization:

(      var_name                         : key name for returned dictionary
       short option character           : same as getopt
       long option string               : same as getopt
       value type                       : str, int, float or None for no associated value
       default                          : default value or None for no associated default
       option_list                      : a Tuple of strings that are acceptable input values
                                          for a string-type variable. Other input values will
                                          fail if this is present.
       option_explanation_list          : a Tuple of explanations of acceptable values for
                                          the string-type variable.
       explanation_strings              : the remainder of the elements in the Tuple are
                                          strings that will be displayed after the options
                                          in the usage listing.
)

"""

DEFAULT_OPTIONS = [(None,
                    "h", "help",
                    None, None,
                    None, None,
                    "display command usage")]

STRING_TYPES = ("s", "S", "str", "STR", "string", "STRING")
INT_TYPES = ("i", "I", "int", "INT", "integer", "INTEGER")
FLOAT_TYPES = ("f", "F", "float", "FLOAT")
FLAG_ON = ("1", "on", "ON", "true", "TRUE", "flag on", "FLAG ON")
FLAG_OFF = ("0", "off", "OFF", "false", "FALSE", "flag off", "FLAG OFF")

def usage(option_list):
    """
    Prints a friendly help message based on the option list.
    """
    # lop off program description and arguments description
    program_description = option_list[0] or ""
    argument_description = option_list[1]
    args = ""
    if argument_description:
        args = "arguments"
        
    print("Usage: %s [options] %s" % (sys.argv[0], args))
    if program_description:
        print(program_description)
        print()
    if argument_description:
        print("Arguments:")
        print(argument_description)
        print()
    print("Options:")
    for option in option_list[2:]:
        var_name, sopt, lopt, vtype, default, opt_list, oe_list, explanation = option
        print("\t-%s, --%-15s : %s" % (sopt[0], lopt, explanation))
        if opt_list:
            print("\tPossible Input Values:")
            for i in range(len(opt_list)):
                opt = opt_list[i]
                try:
                    exp = oe_list[i]
                except (IndexError, TypeError):
                    exp = None
                    
                if not exp:
                    exp = ""
                else:
                    exp = ": %s" % exp
                print("\t\t%-10s %s" % (opt, exp))
        if default:
            print("\t\tDefault: %s" % default)
        if opt_list:
            print()

def process_options(option_list, input=None, add_defaults=1):
    """
    Used in place of getopt. Takes an option list in the format provided above
    and (optionally) an input arg list. If this is not provided the
    sys.argv[1:] list is used pulling directly from the command line.
    """
    input = input or sys.argv[1:]
    short_opts = ""
    long_opts = []
    opt_dict = {}
    return_dict = {}

    if add_defaults:
        new_opts = option_list[0:2]
        for list in (DEFAULT_OPTIONS, option_list[2:]):
            for opt in list:
                new_opts.append(opt)
        option_list = new_opts
    
    for option in option_list[2:]:
        var_name, sopt, lopt, vtype, default, opt_list, oe_list, explanation = option
        col = ""
        eql = ""
        if vtype and (vtype not in FLAG_ON and vtype not in FLAG_OFF):
            col = ":"
            if lopt and not lopt.endswith("="):
                eql = "="
            
        short_opts = "%s%s%s" % (short_opts, sopt[0], col)
        long_opts.append("%s%s" % (lopt, eql))

        opt_dict["-%s" % sopt] = option
        opt_dict["--%s" % lopt] = option

        if not var_name:
            var_name = lopt
            
        return_dict[var_name] = default

    try:
        options, args = getopt(input, short_opts, long_opts)
    except GetoptError:
        usage(option_list)
        sys.exit(2)
    
    for o, a in options:
        if o in opt_dict.keys():
            option = opt_dict[o]
            var_name, sopt, lopt, vtype, default, opt_list, oe_list, explanation = option
            if sopt == "h":
                # special case
                usage(option_list)
                sys.exit(0)

            if vtype:
                if vtype in STRING_TYPES:
                    value = a
                    if opt_list and value not in opt_list:
                        print("Option Error: -%s, --%s: unacceptable input value %s" % (sopt, lopt, value))
                        usage(option_list)
                        sys.exit(2)
                elif vtype in FLAG_ON:
                    value = default or 1
                elif vtype in FLAG_OFF:
                    value = default or 0
                else:
                    if vtype in INT_TYPES:
                        conv = int
                        retstring = "an integer."
                    elif vtype in FLOAT_TYPES:
                        conv = float
                        retstring = "a floating point number."
                    try:
                        value = conv(a)
                    except (TypeError, ValueError):
                        print("Option Error: -%s, --%s must be %s" % (sopt, lopt, retstring))
                        usage(option_list)
                        sys.exit(2)
            else:
                # this defaults to a flag type which is turned on
                value = 1

            if not var_name:
                var_name = lopt

            return_dict[var_name] = value
            
        else:
            # we have a weird problem
            print("Options Internal Error: Unrecognized option %s" % o)
            usage(option_list)
            sys.exit(2)
    
    return_dict["_option_list"] = option_list
    if not args:
        args = []
    #args = [args,]
    
    return return_dict, args
    
            
            
    
                
        






