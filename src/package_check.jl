#package_check.jl

#list includes all packages required for this project
list = String["JuMP", "Gurobi", "BangBang", "SparseArrays", "Random"]

#def function to check installations for all packages in list
function check_packages_installation(list::Vector{String})
    #if list is empty, return a warning
    (length(list) == 0) && (@warn "The list is empty.")
    
    #initialize a notification variable to show whether all packages have been installed
    has_installed = true

    #check packages by for loop
    for package_name in list
        #Base.find_package(package_name) return package path if this package has been installed
        #otherwise it will retun nothing
        if isnothing(Base.find_package(package_name))
            @warn "Please install package: " * string(package_name)
            #update the notification variable
            has_installed = false
        end
    end

    #output
    return has_installed
end

begin
    @info "Checking installations of all necessary pacakges for this project."
    #check installation status
    if check_packages_installation(list) == false
        #error info, we comment next line if we will use try catch control flow
        #@error "There are uninstalled necessary packages. Please install them and check again."
        #for furture rethrow error
        throw(error("There are uninstalled necessary packages. Please install them and check again."))
    else
        @info "All necessary Packages have been installed."
    end
end

#=
We can use a control flow to break program if there are uninstalled necessary packages
try
    include("package_check.jl") #here take an attention to the path of file
catch e
    #this error will interpert the program
    rethrow(e)
end
println("All necessary packages have been installed successfully.")
=#
