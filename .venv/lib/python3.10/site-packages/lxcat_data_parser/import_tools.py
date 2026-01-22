import pandas as pd
from enum import IntEnum
import logging
import os
import io
from typing import Union


# The different types of cross sections
class CrossSectionTypes(IntEnum):
    ELASTIC = 0
    EFFECTIVE = 1
    EXCITATION = 2
    ATTACHMENT = 3
    IONIZATION = 4


CST = CrossSectionTypes


class CrossSection:
    """Class containing data of a single cross section.

    Attributes:
        cross_section_type: Type of collision, possible str/member of CrossSectionTypes
        species: name of the species as a str, example: N2
        data: pandas DataFrame with columns "energy" in eV and "cross section" in m2
        mass_ratio: ratio of electron mass to atomic/molecular mass
        threshold: cross section threshold in eV
        info: optional additional information on the cross section given via kwargs"""

    def __init__(self, cross_section_type: Union[str, CST], species: str,
                 data: pd.DataFrame, mass_ratio: float = None,
                 threshold: float = None, **kwargs):
        if isinstance(cross_section_type, CST):
            self.type = cross_section_type
        else:
            try:
                self.type = CST[cross_section_type]
            except KeyError:
                names = ",".join([xsec.name for xsec in CrossSectionTypes])
                logging.error(f"Invalid value of argument cross_section_type. "
                              f" Accepted values are: {names}")
                raise
        self.species = species
        self.data = data
        self.mass_ratio = mass_ratio
        self.threshold = threshold
        self.info = {}
        for key, value in kwargs.items():
            self.info[key] = value

    def __repr__(self):
        if self.threshold is not None:
            name = (f"{self.species} {self.type.name} CrossSection "
                    f"at {self.threshold} eV")
            return name
        else:
            return f"{self.species} {self.type.name} CrossSection"

    def __eq__(self, other):
        if not isinstance(other, CrossSection):
            return NotImplemented
        if self.type != other.type:
            logging.debug(f"Not the same type: {self.type.name} vs {other.type.name}.")
            return False
        if self.species != other.species:
            logging.debug(f"Not the same species: {self.species} vs {other.species}.")
            return False
        if not self.data.equals(other.data):
            logging.debug("Not the same data.")
            return False
        if self.mass_ratio != other.mass_ratio:
            logging.debug(f"Not the same mass ratio: {self.mass_ratio} vs "
                          f"{other.mass_ratio}.")
            return False
        if self.threshold != other.threshold:
            logging.debug(f"Not the same threshold: {self.threshold} vs "
                          f"{other.threshold}.")
            return False
        if self.info != other.info:
            logging.debug(f"Not the same info: {self.info} vs {other.info}.")
            return False
        return True


class CrossSectionSet:
    """A class containing a set of cross sections.

    Attributes:
        species: name of the species as a str, example: N2
        database: name of the database
        cross_sections: list of CrossSections"""

    def __init__(self, input_file=None, imposed_species=None, imposed_database=None):
        """
        Reads a set of cross sections from a file.

        By default, reads the first cross section set found in the input file but, if
        an imposed_species and/or an imposed_database are defined, reads the first cross
        section set of that species and/or that database found in the input file.
        The input file should be compatible with the LXCat cross section data format.
        """

        self.species = imposed_species
        self.database = imposed_database
        self.cross_sections = []

        if input_file is not None:
            with open(input_file, "r") as f:
                logging.info(f"Starting to read the contents "
                             f"of {os.path.basename(input_file)}")
                current_database = None
                line = f.readline()
                while line:
                    # find the name of the database (optional)
                    if line.startswith("DATABASE:"):
                        current_database = line[9:].strip()
                        line = f.readline()
                    # find a line starting with one of the cross_section_types
                    found_cs = [x.name for x in CST if line.startswith(x.name)]
                    if found_cs:
                        # type of cross section
                        cs_type = CST[found_cs[0]]
                        # species (may be followed by other text on the same line)
                        line = f.readline()
                        current_species = line.split()[0]
                        if self.species is None:
                            self.species = current_species
                        # if this is the right species, proceed
                        if current_species == self.species:
                            if self.database is None:
                                self.database = current_database
                            # if this is the right database, proceed
                            if current_database == self.database:
                                # depending on the type of cross section, the next line
                                # contains either the mass_ratio or the threshold
                                mass_ratio = None
                                threshold = None
                                if cs_type == CST.EFFECTIVE or cs_type == CST.ELASTIC:
                                    mass_ratio = float(f.readline().split()[0])
                                elif any([cs_type == CST.EXCITATION,
                                          cs_type == CST.IONIZATION]):
                                    threshold = float(f.readline().split()[0])
                                # the next lines may contain optional, additional
                                # information on the cross section with the format
                                # KEY: information
                                other_info = {}
                                line = f.readline()
                                while not line.startswith("-----"):
                                    s = line.split(":")
                                    key = s[0].strip()
                                    other_info[key] = line[len(key) + 1:].strip()
                                    line = f.readline()
                                # "-----" mars the start of the tabulated data
                                # put the data into an ioString
                                data_stream = io.StringIO()
                                line = f.readline()
                                while not line.startswith("-----"):
                                    data_stream.write(line)
                                    line = f.readline()
                                data_stream.seek(0)
                                # "-----" marks the end of the tabulated data
                                # read the data into a pandas DataFrame
                                data = pd.read_csv(data_stream, sep="\t",
                                                   names=["energy", "cross section"])
                                # create the cross section object with all the info
                                xsec = CrossSection(cs_type, current_species, data,
                                                    mass_ratio, threshold, **other_info)
                                self.cross_sections.append(xsec)
                    line = f.readline()
                if self.cross_sections:
                    logging.info(f"Initialized {self}")
                else:
                    required = " ".join(s for s in [imposed_database, imposed_species]
                                        if s is not None)
                    logging.error(f"Could not find {required} cross sections "
                                  f"in {os.path.basename(input_file)}")
                    raise CrossSectionReadingError
        else:
            logging.info(f"Initialized {self}")

    def __repr__(self):
        if self.database is not None and self.species is not None:
            return f"{self.database} {self.species} CrossSectionSet"
        elif self.database is not None:
            return f"{self.database} CrossSectionSet"
        elif self.species is not None:
            return f"{self.species} CrossSectionSet"
        else:
            return "Empty CrossSectionSet"

    def write(self, output_file):
        """
        Writes the set of cross sections in a "*.txt" file under "input_file",
        in an LXCat-compatible format.
        """
        with open(output_file, "w") as fh:
            fh.write("Data printed using the package 'lxcat_data_parser', formatted "
                     " in accordance with the cross section data format of LXCat, "
                     "www.lxcat.net.\n\n")
            fh.write("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n")
            if self.database is not None:
                fh.write("DATABASE: " + self.database + "\n")
                fh.write("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n\n")
            fh.write("********************************\n")
            for xsec in self.cross_sections:
                fh.write(xsec.type.name + "\n")
                fh.write(xsec.species + "\n")
                if xsec.type in {CST.ELASTIC, CST.EFFECTIVE}:
                    fh.write(str(xsec.mass_ratio) + "\n")
                elif xsec.type in {CST.EXCITATION, CST.IONIZATION}:
                    fh.write(str(xsec.threshold) + "\n")
                for key in xsec.info.keys():
                    fh.write(key + ": " + xsec.info[key] + "\n")
                # create a 2-column table: "energy" and "values"
                fh.write("-----------------------------\n")
                xsec.data.to_csv(fh, sep="\t", index=False, header=False, chunksize=2,
                                 float_format="%.6e", line_terminator="\n")
                fh.write("-----------------------------\n\n")
            fh.write("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n")

    def __eq__(self, other):
        if not isinstance(other, CrossSectionSet):
            return NotImplemented
        if self.database != other.database:
            return False
        if self.species != other.species:
            return False
        # check that the cross sections are identical (order may vary) by removing them
        # one by one and checking that the remaining list is then empty
        other_xsecs = list(other.cross_sections)
        for xsec in self.cross_sections:
            try:
                other_xsecs.remove(xsec)
            except ValueError:
                return False
        return not other_xsecs


class CrossSectionReadingError(Exception):
    """
    Error indicating a problem with the input file content. Check that the file
    follows the LXCat cross section data format (see www.lxcat.net), and that the
    species and/or database names are correct.
    """
