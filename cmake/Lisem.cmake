# Lisem depends on 3rd party software (eg: Qt, Qwt). Such software can be
# installed by downloading installers or built locally. The Peacock project
# aids in building 3rd party software locally for various platforms. To ease
# the use of software built by Peacock, we include one of its modules and
# use variables defined by it. For more information, visit
# https://github.com/geoneric/peacock

# https://github.com/geoneric/peacock/blob/master/cmake/PeacockPlatform.cmake
INCLUDE(PeacockPlatform) # This one first. Other modules use the variables.

INCLUDE(LisemConfiguration)
INCLUDE(LisemCompiler)
INCLUDE(LisemExternal)
INCLUDE(LisemMacro)
