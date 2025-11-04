from qgis.PyQt.QtCore import QCoreApplication
from qgis.core import QgsProcessingAlgorithm, QgsProcessingProvider


class LinkscapeLaunchAlgorithm(QgsProcessingAlgorithm):
    """Processing entry that opens the Linkscape dialog."""

    ALG_ID = "linkscape_corridor_analysis"

    def __init__(self, plugin):
        super().__init__()
        self._plugin = plugin

    def name(self) -> str:
        return self.ALG_ID

    def displayName(self) -> str:
        return self.tr("Linkscape Corridor Analysis")

    def group(self) -> str:
        return ""

    def groupId(self) -> str:
        return ""

    def shortHelpString(self) -> str:
        return self.tr("Opens the Linkscape dialog so you can configure and run analyses.")

    def flags(self):
        return super().flags() | QgsProcessingAlgorithm.FlagNoThreading

    def tr(self, string: str) -> str:
        return QCoreApplication.translate("LinkscapeLaunchAlgorithm", string)

    def initAlgorithm(self, config=None):  # noqa: D401 - required override
        """No parameters to declare."""
        return

    def processAlgorithm(self, parameters, context, feedback):
        feedback.pushInfo(self.tr("Opening Linkscape dialog..."))
        self._plugin.run()
        return {}

    def createInstance(self):
        return LinkscapeLaunchAlgorithm(self._plugin)


class LinkscapeProcessingProvider(QgsProcessingProvider):
    """Processing provider that exposes the Linkscape launcher."""

    def __init__(self, plugin):
        super().__init__()
        self._plugin = plugin

    def id(self) -> str:
        return "linkscape"

    def name(self) -> str:
        return self.tr("SORUS")

    def tr(self, string: str) -> str:
        return QCoreApplication.translate("LinkscapeProcessingProvider", string)

    def loadAlgorithms(self):
        self.addAlgorithm(LinkscapeLaunchAlgorithm(self._plugin))
