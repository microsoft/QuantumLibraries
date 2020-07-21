// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

#nullable enable

using Microsoft.Jupyter.Core;
using Newtonsoft.Json;

namespace Microsoft.Quantum.Diagnostics
{

    internal class FailureRecord<T>
    {
        [JsonProperty("actual")]
        public T Actual { get; set; }

        [JsonProperty("expected")]
        public T Expected { get; set; }

        [JsonProperty("message")]
        public string Message { get; set; } = "";
    }

    
    public class FailureRecordToHtmlEncoder : IResultEncoder
    {
        public string MimeType => MimeTypes.Html;

        public EncodedData? Encode(object displayable) =>
            displayable is FailureRecord<dynamic> record
            ? $@"
                <blockquote style=""border-left-color: rgb(255, 0, 0)"">
                    <p>
                        <strong>Contradiction reached:</strong> {record.Message}
                    </p>

                    <table>
                        <tr>
                            <td>Actual</td>
                            <td>{record.Actual}</td>
                        </tr>
                        <tr>
                            <td>Expected</td>
                            <td>{record.Expected}</td>
                        </tr>
                    </table>
                </blockquote>                
            ".ToEncodedData()
            : (EncodedData?)null;
    }

    public class FailureRecordToTextEncoder : IResultEncoder
    {
        public string MimeType => MimeTypes.PlainText;

        public EncodedData? Encode(object displayable) =>
            displayable is FailureRecord<dynamic> record
            ? $"{record.Message}\n\tExpected:\t{record.Expected}\n\tActual:\t{record.Actual}".ToEncodedData()
            : (EncodedData?)null;
    }


}
