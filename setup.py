from distutils.core import setup

setup(
    name="ecpy",
    version="1.1",
    description="A Elliptic-Curve Library",
    author="@elliptic_shiho",
    author_email="shiho.elliptic@gmail.com",
    url="https://github.com/elliptic-shiho/ecpy/",
    packages=["ecpy", "ecpy.fields", "ecpy.rings", "ecpy.utils", "ecpy.elliptic_curve"],
    install_requires=["six"],
)
