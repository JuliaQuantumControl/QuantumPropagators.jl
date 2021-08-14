.PHONY: help test docs clean distclean devrepl
.DEFAULT_GOAL := help

define PRINT_HELP_PYSCRIPT
import re, sys

for line in sys.stdin:
    match = re.match(r'^([a-z0-9A-Z_-]+):.*?## (.*)$$', line)
    if match:
        target, help = match.groups()
        print("%-20s %s" % (target, help))
print("""
Instead of "make test", consider "make devrepl" if you want to run the test
suite or generate the docs repeatedly.

Make sure you have Revise.jl installed in your standard Julia environment
""")
endef
export PRINT_HELP_PYSCRIPT


help:  ## show this help
	@python -c "$$PRINT_HELP_PYSCRIPT" < $(MAKEFILE_LIST)


test/Manifest.toml: test/Project.toml
	julia --project=test -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'


test:  ## Run the test suite
	@rm -f test/Manifest.toml  # Pkg.test cannot handle existing Manifest.toml
	julia --startup-file=yes -e 'using Pkg;Pkg.activate(".");Pkg.test(coverage=true)'
	@echo "Done. Consider using 'make devrepl'"


devrepl: test/Manifest.toml ## Start an interactive REPL for testing and building documentation
	@julia --project=test --banner=no --startup-file=yes -e 'include("test/init.jl")' -i


docs/Manifest.toml: docs/Project.toml
	julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'


docs: docs/Manifest.toml ## Build the documentation
	julia --project=docs docs/make.jl
	@echo "Done. Consider using 'make devrepl'"


clean: ## Clean up build/doc/testing artifacts
	rm -f src/*.cov test/*.cov
	rm -f test/examples/*
	for file in examples/*.jl; do rm -f docs/src/"$${file%.jl}".*; done
	rm -rf docs/build


distclean: clean ## Restore to a clean checkout state
	rm -f Manifest.toml docs/Manifest.toml test/Manifest.toml
