import os

from hypothesis import Verbosity, settings

settings.register_profile("ci", max_examples=2000)
settings.register_profile("dev", max_examples=200)
settings.register_profile("debug", max_examples=200, verbosity=Verbosity.verbose)
settings.load_profile(os.getenv(u"HYPOTHESIS_PROFILE", "ci"))
