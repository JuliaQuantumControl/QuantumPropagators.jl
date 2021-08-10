.PHONY: help test docs clean
.DEFAULT_GOAL := help

define PRINT_HELP_PYSCRIPT
import re, sys

for line in sys.stdin:
    match = re.match(r'^([a-z0-9A-Z_-]+):.*?## (.*)$$', line)
    if match:
        target, help = match.groups()
        print("%-20s %s" % (target, help))
print("""

Consider using a long-running Julia REPL for repeated tasks.
""")
endef
export PRINT_HELP_PYSCRIPT


help:  ## show this help
	@python -c "$$PRINT_HELP_PYSCRIPT" < $(MAKEFILE_LIST)


test: ## Run the test suite (pkg> activate .; test)
	julia -e 'using Pkg;Pkg.activate(".");Pkg.test(coverage=true)'

docs: ## Build the documentation
	julia --project=docs docs/make.jl

clean: ## Clean up build/doc/testing artifacts
	rm -f test/examples/*
	for file in examples/*.jl; do rm -f docs/src/"$${file%.jl}".*; done
	rm -rf docs/build
