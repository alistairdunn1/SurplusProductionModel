
PKGNAME = `sed -n "s/Package: *\([^ ]*\)/\1/p" DESCRIPTION`
PKGVERS = `sed -n "s/Version: *\([^ ]*\)/\1/p" DESCRIPTION`

all: clean check

build: dependencies
	R CMD build .

check: build
	R CMD check $(PKGNAME)_$(PKGVERS).tar.gz

dependencies: documentation
	Rscript \
	-e 'if (!requireNamespace("remotes")) install.packages("remotes")' \
	-e 'remotes::install_deps(dependencies = TRUE)'
    
documentation:
	Rscript \
	-e 'if (!requireNamespace("roxygen2")) install.packages("roxygen2")' \
	-e 'roxygen2::roxygenize()'

install: build
	R CMD INSTALL $(PKGNAME)_$(PKGVERS).tar.gz

clean:
	@rm -rf $(PKGNAME)_$(PKGVERS).tar.gz $(PKGNAME).Rcheck
	

